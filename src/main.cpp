// Finds signal / crossing(traffic_signals) / signal triples along highway ways.
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
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <set>
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

void usage(const char* prog) {
    std::cerr << "Usage: " << prog << " [-d meters] [-b meters] [-o out.geojson] input.pbf\n"
              << "  -d  max distance signal<->crossing (default 30)\n"
              << "  -b  bbox buffer (default 10)\n"
              << "  -o  output file (default stdout)\n";
}

} // namespace

int main(int argc, char** argv) {
    double max_dist = 30.0;
    double buffer = 10.0;
    std::string out, in;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-d" && i + 1 < argc) max_dist = std::atof(argv[++i]);
        else if (a == "-b" && i + 1 < argc) buffer = std::atof(argv[++i]);
        else if (a == "-o" && i + 1 < argc) out = argv[++i];
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

        std::FILE* f = out.empty() ? stdout : std::fopen(out.c_str(), "w");
        if (!f) { std::perror(out.c_str()); return 1; }
        write_geojson(f, out_hits);
        if (f != stdout) std::fclose(f);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
