#include <fstream>
#include <string>
#include <regex>
#include <map>
// #include <locale>
#include <jon/full/CPPALPv2.h>
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
    mongocxx::uri uri("mongodb://NEUString:kreuzer@129.10.135.170:27017");
    mongocxx::client conn(uri);

    // Get the collections to use
    mongocxx::collection poly_collection = conn["ToricCY"]["POLY"];
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
        // assert(poly_doc["NALLTRIANGS"].isInt());
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
            LatticePolytope poly = newton_poly.polar();

            std::vector<IntMatrix> face3d_mats = poly.faces(3);

            // std::map<std::string, int> face3d_to_id;
            std::map<std::string, int> face2d_to_id;

            // std::vector<int> face3d_ids;
            // std::vector<int> face2d_ids;

            // std::map<vector<int>, int> face2d_ids_mult;

            std::map<std::string, std::vector<int> > face2d_mat_to_face3d_ids;
            // Create 2FACETO3FACE
            Json::Value face3d_id_to_faces2d_ids;

            // std::vector<IntMatrix> face3d_nf_mats;
            // std::vector<std::vector<IntMatrix> > face3d_face2d_fulldmats;
            // std::vector<IntMatrix> uniq_face2d_fulldmats;
            for (IntMatrix face3d_mat : face3d_mats) {
                IntMatrix cone4d_mat = add_col(face3d_mat, origin_col, 0);
                IntMatrix cone4d_nf_mat = LatticePolytope(cone4d_mat).normal_form();

                std::string cone_doc = jmongo::simple_find_one(cone_collection, "NORMALFORM", intmat_to_string(cone4d_nf_mat));
                int cone_id = std::stoi(jmongo::get_field(cone_doc, "CONEID", reader));
                int cone_triang = std::stoi(jmongo::get_field(cone_doc, "FACETNREGTRIANG", reader));

                // face3d_ids.push_back(cone_id);

                sort_intmat_cols(face3d_mat, compare_cols);

                // face3d_to_id[intmat_to_string(face3d_mat)] = cone_id;

                // IntMatrix face3d_nf_mat = remove_col(cone4d_nf_mat, origin_col);
                // face3d_nf_mats.push_back(face3d_nf_mat);

                // Create the LatticePolytope for the facet
                LatticePolytope face3d_poly(face3d_mat);
                // Get the 2-faces of the facet
                std::vector<IntMatrix> face2d_mats = face3d_poly.faces(2, false);

                // std::vector<int> face2d_ids;

                // std::vector<IntMatrix> face2d_fulldmats;

                Json::Value inner(Json::arrayValue);

                for (IntMatrix face2d_mat : face2d_mats) {
                    // Get 2-face in normal form to obtain FACEID
                    IntMatrix cone3d_nf_mat = LatticePolytope(face2d_mat).normal_form();

                    std::string face_doc = jmongo::simple_find_one(face_collection, "NORMALFORM", intmat_to_string(cone3d_nf_mat));
                    int face_id = std::stoi(jmongo::get_field(face_doc, "FACEID", reader));

                    // face2d_ids.push_back(face_id);

                    inner.append(Json::Value(face_id));

                    // Get 2-face in full dimension, so it can be compared among different 3-faces
                    IntMatrix face2d_fulldmat = face3d_poly.pointsfullD(face2d_mat);
                    // Sort 2-face vertex matrix
                    sort_intmat_cols(face2d_fulldmat, compare_cols);

                    face2d_to_id[intmat_to_string(face2d_fulldmat)] = face_id;

                    face2d_mat_to_face3d_ids[intmat_to_string(face2d_fulldmat)].push_back(cone_id);

                    // face2d_fulldmats.push_back(face2d_fulldmat);

                    // // Add unique 2-faces to uniq_face2d_fulldmats
                    // int i = 0;
                    // while (i < uniq_face2d_fulldmats.size()) {
                    //     if (face2d_fulldmat.cols() == uniq_face2d_fulldmats[i].cols() && face2d_fulldmat == uniq_face2d_fulldmats[i]) { break; }
                    //     i++;
                    // }
                    // if (i == uniq_face2d_fulldmats.size()) {
                    //     uniq_face2d_fulldmats.push_back(face2d_fulldmat);
                    // }
                }

                // std::sort(face2d_ids.begin(), face2d_ids.end());
                // face2d_ids_mult[face2d_ids]++:

                // if (face2d_ids_mult[face2d_ids] == 1) {
                std::string key = std::to_string(cone_id);
                if (face3d_id_to_faces2d_ids.isMember(key)) {
                    face3d_id_to_faces2d_ids[key].append(inner);
                } else {
                    Json::Value outer(Json::arrayValue);
                    outer.append(inner);
                    face3d_id_to_faces2d_ids[key] = outer;
                }
                // }

                // face3d_face2d_fulldmats.push_back(face2d_fulldmats);
            }
 
            // Create 3FACETO2FACE
            // std::map<vector<int>, int> face3d_ids_mult;
            Json::Value face2d_id_to_faces3d_ids;
            for (std::map<std::string, int>::iterator it = face2d_to_id.begin(); it != face2d_to_id.end(); it++) {
                std::string key = std::to_string(it->second);
                std::vector<int> face3d_ids = face2d_mat_to_face3d_ids[it->first];
                // std::sort(face3d_ids.begin(), face3d_ids.end());
                // face3d_ids_mult[face3d_ids]++:
                // if (face3d_ids_mult[face3d_ids] == 1) {
                Json::Value inner(Json::arrayValue);
                for (int face3d_id : face3d_ids) {
                    inner.append(Json::Value(face3d_id));
                }
                if (face2d_id_to_faces3d_ids.isMember(key)) {
                    face2d_id_to_faces3d_ids[key].append(inner);
                } else {
                    Json::Value outer(Json::arrayValue);
                    outer.append(inner);
                    face2d_id_to_faces3d_ids[key] = outer;
                }
                // }
            }

            // std::sort(face3d_ids.begin(), face3d_ids.end());
            // face3d_ids.erase(std::unique(face3d_ids.begin(), face3d_ids.end()), face3d_ids.end());
            // Json::Value face3d_ids_json(Json::arrayValue);
            // for (int face3d_id : face3d_ids) { face3d_ids_json.append(Json::Value(face3d_id)); }

            // std::sort(face2d_ids.begin(), face2d_ids.end());
            // face2d_ids.erase(std::unique(face2d_ids.begin(), face2d_ids.end()), face2d_ids.end());
            // Json::Value face2d_ids_json(Json::arrayValue);
            // for (int face2d_id : face2d_ids) { face2d_ids_json.append(Json::Value(face2d_id)); }

            Json::Value index_doc;
            index_doc["POLYID"] = Json::Value(poly_id);

            Json::Value out_doc;
            out_doc["H11"] = Json::Value(h11);
            out_doc["POLYID"] = Json::Value(poly_id);
            if (poly_doc.isMember("NALLTRIANGS")) {
                int poly_triang = poly_doc["NALLTRIANGS"].asInt();
                out_doc["NFSRT"] = Json::Value(poly_triang);
            }
            // out_doc["CONEIDS"] = face3d_ids_json;
            // out_doc["FACEIDS"] = face2d_ids_json;
            out_doc["CONETOFACE"] = face3d_id_to_faces2d_ids;
            out_doc["FACETOCONE"] = face2d_id_to_faces3d_ids;

            std::cout << "set POLY ";
            writer->write(index_doc, &std::cout);
            std::cout << " ";
            writer->write(out_doc, &std::cout);
            std::cout << std::endl;

            // // Remove duplicate 2-faces
            // for (int i = 0; i < all_face2d_mats.size(); i++) {
            //     for (int j = i + 1; j < all_face2d_mats.size(); j++) {
            //         if (all_face2d_mats[i].cols() == all_face2d_mats[j].cols() && all_face2d_mats[i] == all_face2d_mats[j]) {
            //             all_face2d_mats.erase(all_face2d_mats.begin() + j);
            //             j--;
            //         }
            //     }
            // }

            // double avg = 0;
            // std::map<std::string, std::vector<int> > face2d_to_face3ds;
            // for (IntMatrix face2d_fulldmat : uniq_face2d_fulldmats) {
            //     int count = 0;
            //     for (int i = 0; i < face3d_face2d_fulldmats.size(); i++) {
            //         if (face_in_vec(face3d_face2d_fulldmats[i], face2d_fulldmat)) {
            //             count += LatticePolytope(face3d_mats[i]).p_boundary_points(3).size();
            //         }
            //     }
            //     // cone3d_to_face3d_counts[intmat_to_string(face2d_fulldmat)] = count_face3d;
            //     avg += (double)count - LatticePolytope(face2d_fulldmat).npoints();
            //     // std::cout << cone3d_to_face3d_counts[intmat_to_string(face2d_fulldmat)] << std::endl;
            // }

            // avg /= (double)uniq_face2d_fulldmats.size();

            // std::cout << poly_id << ": " << avg << std::endl;
        }
        std::cout << std::endl;
    }
}