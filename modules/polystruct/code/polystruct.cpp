#include <fstream>
#include <string>
#include <regex>
#include <map>
// #include <locale>
// #include <jon/full/CPPALPv3.h>
#include <CPPALP/helpers.h>
#include <CPPALP/LatticePolytope.h>
#include <jon/full/jpy.h>
#include <jon/full/jmongo.h> //includes jsoncpp
#include <jon/full/jstring.h>

// bsoncxx and mongocxx includes
#include <bsoncxx/builder/stream/document.hpp>
#include <bsoncxx/json.hpp>

#include <mongocxx/client.hpp>
#include <mongocxx/instance.hpp>
#include <mongocxx/stdx.hpp>
#include <mongocxx/uri.hpp>

IntMatrix add_col(const IntMatrix& mat, const IntCol& col, const int pos) {
    if (pos == 0) {
        IntMatrix new_mat(mat.rows(), mat.cols() + 1);
        new_mat.col(0) = col;
        for (int i = 1; i < new_mat.cols(); i++) {
            new_mat.col(i) = mat.col(i - 1);
        }
        return new_mat;
    } else if (pos == 1) {
        IntMatrix new_mat = mat;
        new_mat.conservativeResize(new_mat.rows(), new_mat.cols() + 1);
        new_mat.col(new_mat.cols() - 1) = col;
        return new_mat;
    }
}

IntMatrix remove_col(const IntMatrix& mat, const IntCol& col) {
    std::vector<int> col_inds;
    for (int i = 0; i < mat.cols(); i++)
    {
        if (mat.col(i) != col) col_inds.push_back(i);
    }
    IntMatrix new_mat(mat.rows(), col_inds.size());
    for (int i = 0; i < new_mat.cols(); i++)
    {
        for (int j = 0; j < mat.rows(); j++) {
            new_mat(j, i) = mat(j, col_inds[i]);
        }
    }
    return new_mat;
}

bool face_in_vec(const std::vector<IntMatrix>& face_vec, const IntMatrix& face) {
    bool is_in_vec = false;
    for (int i = 0; i < face_vec.size(); i++) {
        if (face_vec[i].cols() == face.cols() && face_vec[i] == face) {
            is_in_vec = true;
            break;
        }
    }
    return is_in_vec;
}

bool compare_cols(const std::vector<int>& col1, const std::vector<int>& col2) {
    assert(col1.size() == col2.size());
    int i = 0;
    while (i < col1.size() && col1[i] == col2[i]) { i++; }
    if (i == col1.size()) { return false; }
    return col1[i] < col2[i];
}

void sort_intmat_cols(IntMatrix& mat, bool (*compare)(const std::vector<int>&, const std::vector<int>&)) {
    std::vector<std::vector<int> > vec2d;
    for (int i = 0; i < mat.cols(); i++) {
        IntCol col = mat.col(i);
        vec2d.push_back(std::vector<int>(col.data(), col.data() + col.size()));
    }
    std::sort(vec2d.begin(), vec2d.end(), compare);
    for (int i = 0; i < vec2d.size(); i++) {
        for (int j = 0; j < vec2d[0].size(); j++) {
            mat(j, i) = vec2d[i][j];
        }
    }
}

int main(int argc, char* argv[]) {
    // Create the MongoDB instance
    mongocxx::instance inst{};

    // Connect to the database
    mongocxx::uri uri("mongodb://raltman:kreuzer@129.10.135.170:27017/SUBCONES?socketTimeoutMS=1200000");
    mongocxx::client conn(uri);

    // Get the collections to use
    // mongocxx::collection poly_collection = conn["ToricCY"]["POLY"];
    mongocxx::collection cone4d_collection = conn["SUBCONES"]["CONE"];
    mongocxx::collection cone3d_collection = conn["SUBCONES"]["FACE"];

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
        std::vector<std::vector<int> > poly_verts = string_to_vertex_list(nverts);

        IntMatrix poly_mat(poly_verts[0].size(), poly_verts.size());
        for (int i = 0; i < poly_verts.size(); i++) {
            for (int j = 0; j < poly_verts[0].size(); j++) {
                poly_mat(j, i) = (double) poly_verts[i][j];
            }
        }

        IntCol origin_col(poly_verts[0].size());
            for (int i = 0; i < poly_verts[0].size(); i++) {
                origin_col(i) = (double) 0;
            }
        
        LatticePolytope newton_poly = LatticePolytope(poly_mat);

        if (newton_poly.is_reflexive()) {
            // For the dual
            LatticePolytope dual_poly = newton_poly.polar();

            std::vector<IntMatrix> face3d_mats = dual_poly.faces(3);

            std::map<std::string, int> face2d_to_id;
            std::map<std::string, std::vector<int> > face2d_mat_to_face3d_ids;
            // Create 2FACETO3FACE
            std::map<int, int> face3d_id_to_face2d_ids;
            for (IntMatrix face3d_mat : face3d_mats) {
                IntMatrix cone4d_mat = add_col(face3d_mat, origin_col, 0);
                IntMatrix cone4d_nf_mat = LatticePolytope(cone4d_mat).normal_form();

                std::string cone4d_doc = jmongo::simple_find_one(cone4d_collection, "NORMALFORM", intmat_to_string(cone4d_nf_mat));
                int face3d_id = std::stoi(jmongo::get_field(cone4d_doc, "CONEID", reader));

                sort_intmat_cols(face3d_mat, compare_cols);

                // Create the LatticePolytope for the facet
                LatticePolytope face3d_poly(face3d_mat);
                // Get the 2-faces of the facet
                std::vector<IntMatrix> face2d_mats = face3d_poly.faces(2, false);

                // std::vector<int> this_face2d_ids;
                for (IntMatrix face2d_mat : face2d_mats) {
                    // Get 2-face in normal form to obtain FACEID
                    IntMatrix cone3d_nf_mat = LatticePolytope(face2d_mat).normal_form();

                    std::string cone3d_doc = jmongo::simple_find_one(cone3d_collection, "NORMALFORM", intmat_to_string(cone3d_nf_mat));
                    int face2d_id = std::stoi(jmongo::get_field(cone3d_doc, "FACEID", reader));

                    // Get 2-face in full dimension, so it can be compared among different 3-faces
                    IntMatrix face2d_fulldmat = face3d_poly.pointsfullD(face2d_mat);
                    // Sort 2-face vertex matrix
                    sort_intmat_cols(face2d_fulldmat, compare_cols);

                    face2d_to_id[intmat_to_string(face2d_fulldmat)] = face2d_id;

                    face2d_mat_to_face3d_ids[intmat_to_string(face2d_fulldmat)].push_back(face3d_id);
                }

                ++face3d_id_to_face2d_ids[face3d_id];
            }

            std::map<int, std::map<std::string, int> > face2d_id_to_face3d_ids;
            for (std::map<std::string, int>::iterator it1 = face2d_to_id.begin(); it1 != face2d_to_id.end(); it1++) {
                int face2d_id = it1->second;
                std::vector<int> face3d_ids = face2d_mat_to_face3d_ids[it1->first];
                std::vector<int> this_face3d_ids;
                for (int face3d_id : face3d_ids) {
                    this_face3d_ids.push_back(face3d_id);
                }

                std::sort(this_face3d_ids.begin(), this_face3d_ids.end());

                ++face2d_id_to_face3d_ids[face2d_id][vector_to_string(this_face3d_ids)];
            }

            Json::Value json_face3d_id_to_face2d_ids;
            for (std::map<int, int>::iterator it1 = face3d_id_to_face2d_ids.begin(); it1 != face3d_id_to_face2d_ids.end(); ++it1) {
                int face3d_id = it1->first;
                int mult = it1->second;
                std::string key = std::to_string(face3d_id);
                json_face3d_id_to_face2d_ids[key] = Json::Value(mult);
            }

            Json::Value json_face2d_id_to_face3d_ids;
            for (std::map<int, std::map<std::string, int> >::iterator it1 = face2d_id_to_face3d_ids.begin(); it1 != face2d_id_to_face3d_ids.end(); ++it1) {
                int face2d_id = it1->first;
                std::map<std::string, int> face3d_ids_mult = it1->second;
                Json::Value face2d_id_arr(Json::arrayValue);
                for (std::map<std::string, int>::iterator it2 = face3d_ids_mult.begin(); it2 != face3d_ids_mult.end(); ++it2) {
                    std::vector<int> face3d_ids = string_to_vector(it2->first);
                    int mult = it2->second;
                    Json::Value json_face3d_ids_mult(Json::arrayValue);
                    for (int face3d_id : face3d_ids) {
                        json_face3d_ids_mult.append(Json::Value(face3d_id));
                    }
                    json_face3d_ids_mult.append(Json::Value(mult));
                    face2d_id_arr.append(json_face3d_ids_mult);
                }
                std::string key = std::to_string(face2d_id);
                json_face2d_id_to_face3d_ids[key] = face2d_id_arr;
            }

            Json::Value index_doc;
            index_doc["POLYID"] = Json::Value(poly_id);

            Json::Value out_doc;
            out_doc["H11"] = Json::Value(h11);
            out_doc["POLYID"] = Json::Value(poly_id);
            out_doc["NVERTS"] = Json::Value(nverts);

            if (poly_doc.isMember("NALLTRIANGS")) {
                int poly_triang = poly_doc["NALLTRIANGS"].asInt();
                out_doc["NFSRT"] = Json::Value(poly_triang);
            }

            out_doc["CONETOFACE"] = json_face3d_id_to_face2d_ids;
            out_doc["FACETOCONE"] = json_face2d_id_to_face3d_ids;

            std::cout << "set POLY ";
            writer->write(index_doc, &std::cout);
            std::cout << " ";
            writer->write(out_doc, &std::cout);
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}