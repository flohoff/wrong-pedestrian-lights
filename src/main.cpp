// Finds signal / crossing(traffic_signals) / signal triples along highway ways.
#include <sqlite3.h>
#include <spatialite.h>

#include <osmium/geom/haversine.hpp>
#include <osmium/handler.hpp>
#include <osmium/io/pbf_input.hpp>
#include <osmium/io/reader.hpp>
#include <osmium/osm/location.hpp>
#include <osmium/osm/node.hpp>
#include <osmium/osm/way.hpp>
#include <osmium/visitor.hpp>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace {

enum Kind : uint8_t { Signal = 1, Crossing = 2, OtherCrossing = 3 };

struct Poi {
    osmium::Location loc;
    Kind kind;
};

using PoiMap = std::unordered_map<osmium::object_id_type, Poi>;

// Pass 1: collect tagged nodes.
struct NodeHandler : osmium::handler::Handler {
    PoiMap pois;

    void node(const osmium::Node& n) {
        const char* hw = n.tags()["highway"];
        if (!hw || !n.location().valid()) return;
        if (!std::strcmp(hw, "traffic_signals")) {
            pois[n.id()] = {n.location(), Signal};
        } else if (!std::strcmp(hw, "crossing")) {
            const char* c = n.tags()["crossing"];
            bool sig = c && !std::strcmp(c, "traffic_signals");
            pois[n.id()] = {n.location(), sig ? Crossing : OtherCrossing};
        }
    }
};

struct Hit {
    osmium::object_id_type way, a, c, b;
    double min_lon, min_lat, max_lon, max_lat;
    double dist_ac, dist_cb;
    std::vector<osmium::object_id_type> chain; // node ids from a to b
};

bool car_road(const char* hw) {
    static const std::set<std::string> ok = {
        "motorway", "trunk", "primary", "secondary", "tertiary", "unclassified",
        "residential", "service", "living_street", "road", "motorway_link",
        "trunk_link", "primary_link", "secondary_link", "tertiary_link"};
    return ok.count(hw) > 0;
}

struct Item { osmium::object_id_type id; const Poi* poi; size_t pos; };

// Pass 2: scan car ways for signal / one crossing / signal candidates.
struct WayHandler : osmium::handler::Handler {
    const PoiMap& pois;
    double max_dist;
    double buffer;
    std::vector<Hit> hits;
    std::set<std::tuple<osmium::object_id_type, osmium::object_id_type, osmium::object_id_type>> seen;

    WayHandler(const PoiMap& p, double d, double b) : pois(p), max_dist(d), buffer(b) {}

    void way(const osmium::Way& w) {
        const char* hw = w.tags()["highway"];
        if (!hw || !car_road(hw)) return;

        std::vector<Item> items;
        size_t pos = 0;
        for (const auto& nr : w.nodes()) {
            auto it = pois.find(nr.ref());
            if (it != pois.end()) items.push_back({nr.ref(), &it->second, pos});
            ++pos;
        }

        for (size_t i = 0; i < items.size(); ++i) {
            if (items[i].poi->kind != Crossing) continue;
            size_t p = i, n = i;
            bool fp = false, fn = false;
            for (size_t j = i; j-- > 0;)
                if (items[j].poi->kind == Signal) { p = j; fp = true; break; }
            for (size_t j = i + 1; j < items.size(); ++j)
                if (items[j].poi->kind == Signal) { n = j; fn = true; break; }
            if (!fp || !fn || items[p].id == items[n].id) continue;

            // exactly one highway=crossing between the signals
            int crossings = 0;
            for (size_t j = p + 1; j < n; ++j)
                if (items[j].poi->kind != Signal) ++crossings;
            if (crossings != 1) continue;

            const auto& lc = items[i].poi->loc;
            double d1 = osmium::geom::haversine::distance(items[p].poi->loc, lc);
            double d2 = osmium::geom::haversine::distance(lc, items[n].poi->loc);
            if (d1 > max_dist || d2 > max_dist) continue;

            if (!seen.insert({items[p].id, items[i].id, items[n].id}).second) continue;

            const osmium::Location* ls[3] = {&items[p].poi->loc, &lc, &items[n].poi->loc};
            Hit h{w.id(), items[p].id, items[i].id, items[n].id, 1e9, 1e9, -1e9, -1e9, d1, d2, {}};
            for (auto* l : ls) {
                h.min_lon = std::min(h.min_lon, l->lon());
                h.max_lon = std::max(h.max_lon, l->lon());
                h.min_lat = std::min(h.min_lat, l->lat());
                h.max_lat = std::max(h.max_lat, l->lat());
            }
            const double dlat = buffer / 111320.0;
            const double dlon = buffer / (111320.0 * std::cos((h.min_lat + h.max_lat) / 2 * M_PI / 180));
            h.min_lon -= dlon; h.max_lon += dlon;
            h.min_lat -= dlat; h.max_lat += dlat;
            for (size_t k = items[p].pos; k <= items[n].pos; ++k) h.chain.push_back(w.nodes()[k].ref());
            hits.push_back(std::move(h));
        }
    }
};

// Pass 3: neighbours of chain nodes on car ways; >2 distinct neighbours = junction.
struct JunctionHandler : osmium::handler::Handler {
    std::unordered_map<osmium::object_id_type, std::vector<osmium::object_id_type>> nbr;

    explicit JunctionHandler(const std::vector<Hit>& hits) {
        for (const auto& h : hits)
            for (auto id : h.chain) nbr[id];
    }

    void add(osmium::object_id_type id, osmium::object_id_type other) {
        auto it = nbr.find(id);
        if (it == nbr.end()) return;
        auto& v = it->second;
        if (std::find(v.begin(), v.end(), other) == v.end()) v.push_back(other);
    }

    void way(const osmium::Way& w) {
        const char* hw = w.tags()["highway"];
        if (!hw || !car_road(hw)) return;
        const auto& ns = w.nodes();
        for (size_t i = 0; i < ns.size(); ++i) {
            if (i > 0) add(ns[i].ref(), ns[i - 1].ref());
            if (i + 1 < ns.size()) add(ns[i].ref(), ns[i + 1].ref());
        }
    }

    bool junction(osmium::object_id_type id) const {
        auto it = nbr.find(id);
        return it != nbr.end() && it->second.size() > 2;
    }
};

void write_geojson(std::FILE* f, const std::vector<Hit>& hits) {
    std::fputs("{\"type\":\"FeatureCollection\",\"features\":[\n", f);
    for (size_t i = 0; i < hits.size(); ++i) {
        const Hit& h = hits[i];
        std::fprintf(f,
            "{\"type\":\"Feature\",\"properties\":{\"way\":%lld,\"signal_a\":%lld,"
            "\"crossing\":%lld,\"signal_b\":%lld,\"dist_a_m\":%.1f,\"dist_b_m\":%.1f},"
            "\"geometry\":{\"type\":\"Polygon\",\"coordinates\":[[[%.7f,%.7f],[%.7f,%.7f],"
            "[%.7f,%.7f],[%.7f,%.7f],[%.7f,%.7f]]]}}%s\n",
            (long long)h.way, (long long)h.a, (long long)h.c, (long long)h.b, h.dist_ac, h.dist_cb,
            h.min_lon, h.min_lat, h.max_lon, h.min_lat, h.max_lon, h.max_lat,
            h.min_lon, h.max_lat, h.min_lon, h.min_lat, i + 1 < hits.size() ? "," : "");
    }
    std::fputs("]}\n", f);
}

bool valid_ident(const std::string& s) {
    if (s.empty() || std::isdigit((unsigned char)s[0])) return false;
    return std::all_of(s.begin(), s.end(), [](unsigned char c) { return std::isalnum(c) || c == 0x5f; });
}

void check(int rc, sqlite3* db, const char* what) {
    if (rc != SQLITE_OK && rc != SQLITE_DONE && rc != SQLITE_ROW)
        throw std::runtime_error(std::string(what) + ": " + sqlite3_errmsg(db));
}

void exec(sqlite3* db, const std::string& sql) {
    char* err = nullptr;
    if (sqlite3_exec(db, sql.c_str(), nullptr, nullptr, &err) != SQLITE_OK) {
        std::string m = err ? err : "unknown";
        sqlite3_free(err);
        throw std::runtime_error(m + " [" + sql + "]");
    }
}

void write_spatialite(const std::string& path, const std::string& table, const std::string& geom,
                      const std::vector<Hit>& hits) {
    if (!valid_ident(table) || !valid_ident(geom))
        throw std::runtime_error("table/column names must match [A-Za-z_][A-Za-z0-9_]*");

    sqlite3* db = nullptr;
    if (sqlite3_open(path.c_str(), &db) != SQLITE_OK) {
        std::string m = db ? sqlite3_errmsg(db) : "open failed";
        sqlite3_close(db);
        throw std::runtime_error(m);
    }
    void* cache = spatialite_alloc_connection();
    spatialite_init_ex(db, cache, 0);
    try {
        // fresh file: create spatial metadata
        sqlite3_stmt* st = nullptr;
        check(sqlite3_prepare_v2(db, "SELECT 1 FROM sqlite_master WHERE name=\x27geometry_columns\x27", -1, &st, nullptr), db, "prepare");
        bool has_meta = sqlite3_step(st) == SQLITE_ROW;
        sqlite3_finalize(st);
        if (!has_meta) exec(db, "SELECT InitSpatialMetaData(1)");

        check(sqlite3_prepare_v2(db, "SELECT 1 FROM sqlite_master WHERE name=?", -1, &st, nullptr), db, "prepare");
        sqlite3_bind_text(st, 1, table.c_str(), -1, SQLITE_STATIC);
        bool exists = sqlite3_step(st) == SQLITE_ROW;
        sqlite3_finalize(st);
        if (exists) throw std::runtime_error("table \x27" + table + "\x27 already exists in " + path);

        exec(db, "BEGIN");
        exec(db, "CREATE TABLE \"" + table + "\" (id INTEGER PRIMARY KEY AUTOINCREMENT, way INTEGER, "
                 "signal_a INTEGER, crossing INTEGER, signal_b INTEGER, dist_a_m REAL, dist_b_m REAL)");
        exec(db, "SELECT AddGeometryColumn(\x27" + table + "\x27, \x27" + geom + "\x27, 4326, \x27POLYGON\x27, \x27XY\x27)");

        std::string sql = "INSERT INTO \"" + table + "\" (way, signal_a, crossing, signal_b, dist_a_m, dist_b_m, \"" +
                          geom + "\") VALUES (?,?,?,?,?,?,GeomFromText(?,4326))";
        check(sqlite3_prepare_v2(db, sql.c_str(), -1, &st, nullptr), db, "prepare insert");
        for (const Hit& h : hits) {
            char wkt[256];
            std::snprintf(wkt, sizeof wkt, "POLYGON((%.7f %.7f,%.7f %.7f,%.7f %.7f,%.7f %.7f,%.7f %.7f))",
                          h.min_lon, h.min_lat, h.max_lon, h.min_lat, h.max_lon, h.max_lat,
                          h.min_lon, h.max_lat, h.min_lon, h.min_lat);
            sqlite3_bind_int64(st, 1, h.way);
            sqlite3_bind_int64(st, 2, h.a);
            sqlite3_bind_int64(st, 3, h.c);
            sqlite3_bind_int64(st, 4, h.b);
            sqlite3_bind_double(st, 5, h.dist_ac);
            sqlite3_bind_double(st, 6, h.dist_cb);
            sqlite3_bind_text(st, 7, wkt, -1, SQLITE_TRANSIENT);
            check(sqlite3_step(st), db, "insert");
            sqlite3_reset(st);
        }
        sqlite3_finalize(st);
        exec(db, "COMMIT");
    } catch (...) {
        sqlite3_close(db);
        spatialite_cleanup_ex(cache);
        throw;
    }
    sqlite3_close(db);
    spatialite_cleanup_ex(cache);
}

void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " [-d meters] [-b meters] [-o out.geojson] input.pbf\n"
              << "  -d  max distance signal<->crossing (default 30)\n"
              << "  -b  bbox buffer (default 10)\n"
              << "  -o  GeoJSON output file (default stdout unless -s given)\n"
              << "  -s  SpatiaLite output file (fails if table exists)\n"
              << "  -t  SpatiaLite table name (default crossings)\n"
              << "  -g  SpatiaLite geometry column (default geom)\n";
}

} // namespace

int main(int argc, char** argv) {
    double max_dist = 30.0;
    double buffer = 10.0;
    std::string out, in;
    std::string sqlite_out, table = "crossings", geom_col = "geom";
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-d" && i + 1 < argc) max_dist = std::atof(argv[++i]);
        else if (a == "-b" && i + 1 < argc) buffer = std::atof(argv[++i]);
        else if (a == "-o" && i + 1 < argc) out = argv[++i];
        else if (a == "-s" && i + 1 < argc) sqlite_out = argv[++i];
        else if (a == "-t" && i + 1 < argc) table = argv[++i];
        else if (a == "-g" && i + 1 < argc) geom_col = argv[++i];
        else if (a == "-h" || a == "--help") { usage(argv[0]); return 0; }
        else if (in.empty()) in = a;
        else { usage(argv[0]); return 1; }
    }
    if (in.empty()) { usage(argv[0]); return 1; }

    try {
        NodeHandler nh;
        {
            osmium::io::Reader r{in, osmium::osm_entity_bits::node};
            osmium::apply(r, nh);
            r.close();
        }
        std::cerr << "tagged nodes: " << nh.pois.size() << "\n";

        WayHandler wh{nh.pois, max_dist, buffer};
        {
            osmium::io::Reader r{in, osmium::osm_entity_bits::way};
            osmium::apply(r, wh);
            r.close();
        }

        JunctionHandler jh{wh.hits};
        {
            osmium::io::Reader r{in, osmium::osm_entity_bits::way};
            osmium::apply(r, jh);
            r.close();
        }
        std::vector<Hit> out_hits;
        for (auto& h : wh.hits)
            if (std::none_of(h.chain.begin(), h.chain.end(), [&](auto id) { return jh.junction(id); }))
                out_hits.push_back(std::move(h));
        std::cerr << "matches: " << out_hits.size() << "\n";

        if (!sqlite_out.empty()) write_spatialite(sqlite_out, table, geom_col, out_hits);
        if (!out.empty() || sqlite_out.empty()) {
            std::FILE* f = out.empty() ? stdout : std::fopen(out.c_str(), "w");
            if (!f) { std::perror(out.c_str()); return 1; }
            write_geojson(f, out_hits);
            if (f != stdout) std::fclose(f);
        }
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
