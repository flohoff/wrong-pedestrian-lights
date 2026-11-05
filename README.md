# wpl

Find traffic_signals / crossing=traffic_signals / traffic_signals triples on one highway.

    cmake -S . -B build && cmake --build build
    ./build/wpl [-d meters] [-b meters] [-o out.geojson] [-s out.sqlite -t table -g geomcol] in.pbf
