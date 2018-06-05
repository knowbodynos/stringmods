#include <fstream>
#include <string>
#include <regex>
#include <map>
#include <numeric>
#include <cmath>
#include <jon/full/CPPALPv2.h>
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
    // Create the MongoDB instance
    mongocxx::instance inst{};

    // Connect to the database
    mongocxx::uri uri("mongodb://raltman:kreuzer@129.10.135.170:27017/SUBCONES");
    mongocxx::client conn(uri);

    // Get the collections to use
    // mongocxx::collection poly_collection = conn["ToricCY"]["POLY"];
    mongocxx::collection cone_collection = conn["SUBCONES"]["CONE"];
    mongocxx::collection face_collection = conn["SUBCONES"]["FACE"];

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
        assert(poly_doc["NFSRT"].isInt());
        assert(poly_doc["CONETOFACE"].isObject());
        assert(poly_doc["FACETOCONE"].isObject());

        int h11 = poly_doc["H11"].asInt();
        int poly_id = poly_doc["POLYID"].asInt();
        std::string nverts = poly_doc["NVERTS"].asString();
        int nfsrt = poly_doc["NFSRT"].asInt();
        Json::Value conetoface = poly_doc["CONETOFACE"];
        Json::Value facetocone = poly_doc["FACETOCONE"];

        LatticePolytope dual_poly = LatticePolytope(string_to_intmat(nverts)).polar();

        int nfsrt_predict = 1;

        std::map<int, std::map<std::string, double> > cone_info;
        for(Json::Value::iterator it = conetoface.begin(); it != conetoface.end(); ++it) {
            int cone_id = std::stoi(it.key().asString());
            Json::Value cone_face_arr = (*it);
            std::string cone_doc = jmongo::simple_find_one(cone_collection, "CONEID", cone_id);
            cone_info[cone_id]["3FACEN3VERTS"] = std::stod(jmongo::get_field(cone_doc, "NFACETVERTS", reader));
            cone_info[cone_id]["3FACEN3POINTS"] = std::stod(jmongo::get_field(cone_doc, "NFACETPOINTS", reader));
            cone_info[cone_id]["3FACEN2POINTS"] = std::stod(jmongo::get_field(cone_doc, "NSKEL2POINTS", reader));
            cone_info[cone_id]["3FACEN1POINTS"] = std::stod(jmongo::get_field(cone_doc, "NSKEL1POINTS", reader));
            cone_info[cone_id]["3FACE3VOLUME"] = std::stod(jmongo::get_field(cone_doc, "FACETVOLUME", reader));
            cone_info[cone_id]["3FACEN2FACES"] = std::stod(jmongo::get_field(cone_doc, "NFACETFACES", reader));
            cone_info[cone_id]["3FACEN1FACES"] = std::stod(jmongo::get_field(cone_doc, "NFACETEDGES", reader));
            int nfrt = std::stoi(jmongo::get_field(cone_doc, "FACETNREGTRIANG", reader));
            nfsrt_predict *= std::pow(nfrt, cone_face_arr.size());
        }

        std::map<int, std::map<std::string, double> > face_info;
        std::map<int, std::map<std::string, std::vector<double> > > face_cone_info;
        for(Json::Value::iterator it1 = facetocone.begin(); it1 != facetocone.end(); ++it1) {
            int face_id = std::stoi(it1.key().asString());
            std::string face_doc = jmongo::simple_find_one(face_collection, "FACEID", face_id);
            face_info[face_id]["2FACEN2VERTS"] = std::stod(jmongo::get_field(face_doc, "NVERTS", reader));
            face_info[face_id]["2FACEN2POINTS"] = std::stod(jmongo::get_field(face_doc, "NPOINTS", reader));
            face_info[face_id]["2FACEN1POINTS"] = std::stod(jmongo::get_field(face_doc, "NBDPOINTS", reader));
            std::string normal_form = jmongo::get_field(face_doc, "NORMALFORM", reader);
            face_info[face_id]["2FACE2VOLUME"] = LatticePolytope(string_to_intmat(normal_form)).volume();
            face_info[face_id]["2FACEN1FACES"] = std::stod(jmongo::get_field(face_doc, "NEDGES", reader));
            Json::Value face_cone_arr = (*it1);
            for (Json::Value cone_arr : face_cone_arr) {
                std::map<std::string, double> this_face_cone_info;
                for (Json::Value cone_id_json : cone_arr) {
                    int cone_id = cone_id_json.asInt();
                    std::map<std::string, double> this_cone_info = cone_info[cone_id];
                    for(std::map<std::string, double>::iterator it2 = this_cone_info.begin(); it2 != this_cone_info.end(); ++it2) {
                        std::string key = it2->first;
                        double value = it2->second;
                        this_face_cone_info[key] += value;
                    }
                }
                this_face_cone_info["3FACEN3VERTS"] -= face_info[face_id]["2FACEN2VERTS"];
                this_face_cone_info["3FACEN3POINTS"] -= face_info[face_id]["2FACEN2POINTS"];
                this_face_cone_info["3FACEN2POINTS"] -= face_info[face_id]["2FACEN2POINTS"];
                this_face_cone_info["3FACEN1POINTS"] -= face_info[face_id]["2FACEN1POINTS"];
                this_face_cone_info["3FACEN2FACES"] -= 1;
                this_face_cone_info["3FACEN1FACES"] -= face_info[face_id]["2FACEN1FACES"];

                face_cone_info[face_id]["COMPLEXN3VERTS"].push_back(this_face_cone_info["3FACEN3VERTS"]);
                face_cone_info[face_id]["COMPLEXN3POINTS"].push_back(this_face_cone_info["3FACEN3POINTS"]);
                face_cone_info[face_id]["COMPLEXN2POINTS"].push_back(this_face_cone_info["3FACEN2POINTS"]);
                face_cone_info[face_id]["COMPLEXN1POINTS"].push_back(this_face_cone_info["3FACEN1POINTS"]);
                face_cone_info[face_id]["COMPLEX3VOLUME"].push_back(this_face_cone_info["3FACE3VOLUME"]);
                face_cone_info[face_id]["COMPLEXN2FACES"].push_back(this_face_cone_info["3FACEN2FACES"]);
                face_cone_info[face_id]["COMPLEXN1FACES"].push_back(this_face_cone_info["3FACEN1FACES"]);
            }
        }

        // Json::Value cone_avgs;
        std::map<std::string, std::vector<double> > cone_fields;
        for(Json::Value::iterator it1 = conetoface.begin(); it1 != conetoface.end(); ++it1) {
            int cone_id = std::stoi(it1.key().asString());
            Json::Value cone_face_arr = (*it1);
            std::map<std::string, double> this_cone_info = cone_info[cone_id];
            // for (Json::Value face_arr : cone_face_arr) {
            for(std::map<std::string, double>::iterator it2 = this_cone_info.begin(); it2 != this_cone_info.end(); ++it2) {
                std::string key = it2->first;
                double value = it2->second;
                // double new_value = cone_avgs[key].asDouble() + (value * (double)cone_face_arr.size());
                // cone_avgs[key] = Json::Value(new_value);
                for (int i = 0; i < cone_face_arr.size(); i++) {
                    cone_fields[key].push_back(value);
                }
            }
            // }
        }

        // Json::Value face_avgs;
        std::map<std::string, std::vector<double> > face_fields;
        for(Json::Value::iterator it1 = facetocone.begin(); it1 != facetocone.end(); ++it1) {
            int face_id = std::stoi(it1.key().asString());
            Json::Value face_cone_arr = (*it1);
            std::map<std::string, double> this_face_info = face_info[face_id];
            for(std::map<std::string, double>::iterator it2 = this_face_info.begin(); it2 != this_face_info.end(); ++it2) {
                std::string key = it2->first;
                double value = it2->second;
                // double new_value = face_avgs[key].asDouble() + (value * (double)face_cone_arr.size());
                // face_avgs[key] = Json::Value(new_value);
                for (int i = 0; i < face_cone_arr.size(); i++) {
                    face_fields[key].push_back(value);
                }
            }
            std::map<std::string, std::vector<double> > this_face_cone_info = face_cone_info[face_id];
            for(std::map<std::string, std::vector<double> >::iterator it2 = this_face_cone_info.begin(); it2 != this_face_cone_info.end(); ++it2) {
                std::string key = it2->first;
                std::vector<double> values = it2->second;
                for (double value : values) {
                    // double new_value = face_avgs[key].asDouble() + value;
                    // face_avgs[key] = Json::Value(new_value);
                    face_fields[key].push_back(value);
                }
            }
        }
        
        Json::Value index_doc;
        index_doc["POLYID"] = Json::Value(poly_id);

        Json::Value out_doc;
        out_doc["H11"] = Json::Value(h11);
        out_doc["POLYID"] = Json::Value(poly_id);
        out_doc["4POLYN4VERTS"] = dual_poly.nvertices();
        out_doc["4POLYN4POINTS"] = dual_poly.npoints();
        out_doc["4POLYN3POINTS"] = dual_poly.p_boundary_indices(4).size();
        out_doc["4POLYN2POINTS"] = dual_poly.p_boundary_indices(3).size();
        out_doc["4POLYN1POINTS"] = dual_poly.p_boundary_indices(2).size();
        out_doc["4POLYN3FACES"] = dual_poly.faces(3).size();
        out_doc["4POLYN2FACES"] = dual_poly.faces(2).size();
        out_doc["4POLYN1FACES"] = dual_poly.faces(1).size();
        out_doc["NFSRT"] = Json::Value(nfsrt);
        out_doc["NFSRTPREDICT"] = Json::Value(nfsrt_predict);
        for(std::map<std::string, std::vector<double> >::iterator it = cone_fields.begin(); it != cone_fields.end(); ++it) {
            std::string key = it->first;
            std::vector<double> values = it->second;
            double min = *std::min_element(values.begin(), values.end());
            double max = *std::max_element(values.begin(), values.end());
            double sum = std::accumulate(values.begin(), values.end(), 0.0);
            double mean = sum / values.size();
            double sq_sum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
            double stdev = std::sqrt(sq_sum / values.size() - mean * mean);
            out_doc[key]["MIN"] = min;
            out_doc[key]["MAX"] = max;
            out_doc[key]["MEAN"] = mean;
            out_doc[key]["STDEV"] = stdev;
        }
        for(std::map<std::string, std::vector<double> >::iterator it = face_fields.begin(); it != face_fields.end(); ++it) {
            std::string key = it->first;
            std::vector<double> values = it->second;
            double min = *std::min_element(values.begin(), values.end());
            double max = *std::max_element(values.begin(), values.end());
            double sum = std::accumulate(values.begin(), values.end(), 0.0);
            double mean = sum / values.size();
            double sq_sum = std::inner_product(values.begin(), values.end(), values.begin(), 0.0);
            double stdev = std::sqrt(sq_sum / values.size() - mean * mean);
            out_doc[key]["MIN"] = min;
            out_doc[key]["MAX"] = max;
            out_doc[key]["MEAN"] = mean;
            out_doc[key]["STDEV"] = stdev;
        }

        std::cout << "set OVERCOUNT ";
        writer->write(index_doc, &std::cout);
        std::cout << " ";
        writer->write(out_doc, &std::cout);
        std::cout << std::endl << std::endl;
    }
}