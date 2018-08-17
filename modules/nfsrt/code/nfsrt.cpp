#include <iostream>
// #include <fstream>
#include <string>
#include <regex>
#include <map>
// #include <numeric>
// #include <cmath>

// #include <jon/full/CPPALPv3.h>
#include <CPPALP/helpers.h>
#include <CPPALP/LatticePolytope.h>
// #include <jon/full/jpy.h>
#include <jon/full/jmongo.h> //includes jsoncpp
// #include <jon/full/jstring.h>

// bsoncxx and mongocxx includes
#include <bsoncxx/builder/stream/document.hpp>
#include <bsoncxx/json.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>
#include <mongocxx/stdx.hpp>
#include <mongocxx/uri.hpp>

int main(int argc, char* argv[]) {
    // std::string nverts = "{{1,0,0,0},{1,2,0,0},{0,0,1,0},{0,3,1,0},{0,0,0,1},{0,3,0,1},{-1,0,2,-2},{-1,-2,-2,2},{-1,0,-2,2},{1,-2,-2,0},{1,-2,0,-2},{-1,-2,2,-2},{2,-1,-1,-1}}";
    // LatticePolytope poly = LatticePolytope(string_to_intmat(nverts)).polar();

    // std::cout << poly.FSRTs() << std::endl;

    // // Create the MongoDB instance
    // mongocxx::instance inst{};

    // // Connect to the database
    // mongocxx::uri uri("mongodb://raltman:kreuzer@129.10.135.170:27017/SUBCONES");
    // mongocxx::client conn(uri);

    // // Get the collections to use
    // // mongocxx::collection poly_collection = conn["ToricCY"]["POLY"];
    // mongocxx::collection cone_collection = conn["SUBCONES"]["CONE"];
    // mongocxx::collection face_collection = conn["SUBCONES"]["FACE"];

    // Query for the data
    Json::CharReaderBuilder readerBuilder;
    Json::CharReader* reader = readerBuilder.newCharReader();
    Json::StreamWriterBuilder writerBuilder;
    writerBuilder["commentStyle"] = "None";
    writerBuilder["indentation"] = "";
    writerBuilder["enableYAMLCompatibility"] = false;
    std::unique_ptr<Json::StreamWriter> writer(writerBuilder.newStreamWriter());

    Json::Value poly_doc;
    std::string err;
    bool success;

    // std::locale::global(std::locale("en_US.utf8"));

    // std::wstring wline;
    // while (getline(std::wcin, wline)) {
    std::string line;
    while (getline(std::cin, line)) {
    // getline(std::cin, line);
    // By h11
    // int h11a = std::stoi(argv[1]);
    // auto poly_curs = jmongo::simple_find(poly_collection, "H11", h11a);
    // for (auto&& poly_doc_bson : poly_curs) {
    //     // Parse line to json
    //     std::string line = bsoncxx::to_json(poly_doc_bson);
        // std::string line = std::string(wline.begin(), wline.end());
        line = std::regex_replace(line, std::regex("'"), "\"");
        success = reader->parse(line.c_str(), line.c_str() + line.size(), &poly_doc, &err);
        if (!success) {
            std::cout << line << std::endl;
            std::cout << err << std::endl;
        }

        assert(poly_doc["H11"].isInt());
        assert(poly_doc["POLYID"].isInt());
        assert(poly_doc["NVERTS"].isString());

        int h11 = poly_doc["H11"].asInt();
        int poly_id = poly_doc["POLYID"].asInt();
        std::string nverts = poly_doc["NVERTS"].asString();

        IntMatrix poly_mat = string_to_intmat(nverts);

        LatticePolytope dual_poly = LatticePolytope(poly_mat).polar();
        int nfsrt = dual_poly.FSRTs().size();   

        Json::Value index_doc;
        index_doc["POLYID"] = Json::Value(poly_id);

        Json::Value out_doc;

        out_doc["H11"] = Json::Value(h11);
        out_doc["POLYID"] = Json::Value(poly_id);
        out_doc["NVERTS"] = Json::Value(nverts);
        out_doc["NFSRT"] = Json::Value(nfsrt);

        std::cout << "set POLY ";
        writer->write(index_doc, &std::cout);
        std::cout << " ";
        writer->write(out_doc, &std::cout);
        std::cout << std::endl << std::endl;
    }
}