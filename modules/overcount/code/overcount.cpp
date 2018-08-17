#include <iostream>
// #include <fstream>
#include <string>
#include <regex>
#include <map>
// #include <locale>
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
    mongocxx::collection face2d_collection = conn["SUBCONES"]["FACE"];

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

        IntCol origin_col(poly_mat.rows());
        for (int i = 0; i < poly_mat.rows(); i++) {
            origin_col(i) = (double) 0;
        }
        
        LatticePolytope newton_poly = LatticePolytope(poly_mat);

        if (newton_poly.is_reflexive()) {
            // For the dual
            LatticePolytope dual_poly = LatticePolytope(newton_poly.polar().normal_form());

            // std::vector<IntMatrix> test_face2d_mats = dual_poly.faces(2, false);
            // for (IntMatrix test_face2d_mat : test_face2d_mats) {
            //     LatticePolytope test_cone3d_poly = LatticePolytope(test_face2d_mat);
            //     IntMatrix test_face2d_nfmat = test_cone3d_poly.normal_form();
            //     std::cout << intmat_to_string(test_face2d_nfmat) << std::endl;
            // }
            // std::cout << std::endl << std::endl;

            std::vector<IntMatrix> face3d_mats = dual_poly.faces(3, false);

            // std::map<std::string, int> face2d_fulldmat_to_face2d_nfmat;
            std::map<std::string, std::string> face2d_fulldmat_to_face2d_nfmat;

            std::map<std::string, std::vector<std::string> > face2d_fulldmat_to_cone4d_nfmats;
            // Create 2FACETO3FACE
            // Json::Value cone4d_nfmat_to_mult;
            std::map<std::string, int> cone4d_nfmat_to_mult;

            double ln_nfrt_sum = 0;
            double ln_nfrt_predict_sum = 0;

            std::map<std::string, std::map<std::string, double> > face3d_info;
            std::map<std::string, std::string> face3d_nfrt;
            std::map<std::string, double> face3d_nfrt_predict;
            std::map<std::string, std::vector<double> > face3d_fields;
            // std::map<int, std::map<std::string, double> > face2d_info;
            std::map<std::string, std::map<std::string, double> > face2d_info;
            std::map<std::string, std::vector<double> > face2d_fields;
            for (IntMatrix face3d_mat : face3d_mats) {
                // Create the LatticePolytope for the cone and facet
                LatticePolytope cone4d_poly(face3d_mat);
                IntMatrix face3d_nfmat = cone4d_poly.normal_form();
                LatticePolytope face3d_poly(face3d_nfmat);

                IntMatrix cone4d_mat = add_col(face3d_mat, origin_col, 0);
                IntMatrix cone4d_nfmat = LatticePolytope(cone4d_mat).normal_form();

                std::map<std::string, std::map<std::string, double> >::iterator it = face3d_info.find(intmat_to_string(cone4d_nfmat));
                if (it == face3d_info.end()) {
                    // std::cout << "aa: " << intmat_to_string(cone4d_nfmat) << std::endl;
                    std::string cone4d_doc = jmongo::simple_find_one(cone4d_collection, "NORMALFORM", intmat_to_string(cone4d_nfmat));
                    int face3d_id = std::stoi(jmongo::get_field(cone4d_doc, "CONEID", reader));

                    std::map<std::string, double> this_face3d_info;
                    this_face3d_info["3FACEN3VERTS"] = std::stod(jmongo::get_field(cone4d_doc, "NFACETVERTS", reader));
                    this_face3d_info["3FACEN3POINTS"] = std::stod(jmongo::get_field(cone4d_doc, "NFACETPOINTS", reader));
                    this_face3d_info["3FACEN2POINTS"] = std::stod(jmongo::get_field(cone4d_doc, "NSKEL2POINTS", reader));
                    this_face3d_info["3FACEN1POINTS"] = std::stod(jmongo::get_field(cone4d_doc, "NSKEL1POINTS", reader));
                    // this_face3d_info["3FACE3VOLUME"] = std::stod(jmongo::get_field(cone4d_doc, "FACETVOLUME", reader));
                    // std::cout << "a: " << intmat_to_string(cone4d_poly.normal_form()) << " " << face3d_poly.npoints() << std::endl;
                    this_face3d_info["3FACE3VOLUME"] = face3d_poly.volume();
                    this_face3d_info["3FACEN2FACES"] = std::stod(jmongo::get_field(cone4d_doc, "NFACETFACES", reader));
                    this_face3d_info["3FACEN1FACES"] = std::stod(jmongo::get_field(cone4d_doc, "NFACETEDGES", reader));

                    face3d_info[intmat_to_string(cone4d_nfmat)] = this_face3d_info;
                    face3d_nfrt[intmat_to_string(cone4d_nfmat)] = jmongo::get_field(cone4d_doc, "FACETNREGTRIANG", reader);
                    face3d_nfrt_predict[intmat_to_string(cone4d_nfmat)] = std::stod(jmongo::get_field(cone4d_doc, "LNNFRTPREDICT", reader));
                }

                for (std::map<std::string, double>::iterator it = face3d_info[intmat_to_string(cone4d_nfmat)].begin(); it != face3d_info[intmat_to_string(cone4d_nfmat)].end(); ++it) {
                    std::string key = it->first;
                    double value = it->second;
                    // face3d_info[intmat_to_string(cone4d_nfmat)][key] = value;
                    face3d_fields[key].push_back(value);
                }
                
                double ln_nfrt_predict = face3d_nfrt_predict[intmat_to_string(cone4d_nfmat)];
                ln_nfrt_predict_sum += ln_nfrt_predict;

                std::string nfrt_str = face3d_nfrt[intmat_to_string(cone4d_nfmat)];
                if (nfrt_str == "") {
                    ln_nfrt_sum += ln_nfrt_predict;
                } else {
                    int ln_nfrt = std::log(std::stod(nfrt_str));
                    ln_nfrt_sum += ln_nfrt;
                }

                // sort_intmat_cols(face3d_mat, compare_cols);

                // Get the 2-faces of the facet
                std::vector<IntMatrix> face2d_mats = cone4d_poly.faces(2, false);

                // Json::Value inner(Json::arrayValue);

                // std::vector<int> this_face2d_ids;
                for (IntMatrix face2d_mat : face2d_mats) {
                    // Get 2-face in normal form to obtain FACEID
                    LatticePolytope cone3d_poly = LatticePolytope(face2d_mat);
                    IntMatrix face2d_nfmat = cone3d_poly.normal_form();
                    LatticePolytope face2d_poly = LatticePolytope(face2d_nfmat);

                    if (face2d_info.find(intmat_to_string(face2d_nfmat)) == face2d_info.end()) {
                        // std::cout << "bb: " << intmat_to_string(face2d_nfmat) << std::endl;
                        // std::string face2d_doc = jmongo::simple_find_one(face2d_collection, "NORMALFORM", intmat_to_string(face2d_nfmat));
                        // int face2d_id = std::stoi(jmongo::get_field(face2d_doc, "FACEID", reader));

                        std::map<std::string, double> this_face2d_info;
                        // this_face2d_info["2FACEN2VERTS"] = std::stod(jmongo::get_field(face2d_doc, "NVERTS", reader));
                        this_face2d_info["2FACEN2VERTS"] = (double)(face2d_poly.nvertices());
                        // this_face2d_info["2FACEN2POINTS"] = std::stod(jmongo::get_field(face2d_doc, "NPOINTS", reader));
                        this_face2d_info["2FACEN2POINTS"] = (double)(face2d_poly.npoints());
                        // this_face2d_info["2FACEN1POINTS"] = std::stod(jmongo::get_field(face2d_doc, "NBDPOINTS", reader));
                        this_face2d_info["2FACEN1POINTS"] = (double)(face2d_poly.p_boundary_points(2).cols());
                        // std::cout << "b: " << intmat_to_string(face2d_nfmat) << std::endl;
                        this_face2d_info["2FACE2VOLUME"] = face2d_poly.volume();
                        // this_face2d_info["2FACEN1FACES"] = std::stod(jmongo::get_field(face2d_doc, "NEDGES", reader));
                        this_face2d_info["2FACEN1FACES"] = (double)(face2d_poly.faces(1).size());

                        face2d_info[intmat_to_string(face2d_nfmat)] = this_face2d_info;
                    }

                    // Get 2-face in full dimension, so it can be compared among different 3-faces
                    IntMatrix face2d_fulldmat = dual_poly.pointsfullD(face2d_mat);

                    // Sort 2-face vertex matrix
                    sort_intmat_cols(face2d_fulldmat, compare_cols);

                    // std::map<std::string, int>::iterator it = face2d_fulldmat_to_face2d_nfmat.find(intmat_to_string(face2d_fulldmat));
                    std::map<std::string, std::string>::iterator it = face2d_fulldmat_to_face2d_nfmat.find(intmat_to_string(face2d_fulldmat));
                    for (std::map<std::string, double>::iterator it2 = face2d_info[intmat_to_string(face2d_nfmat)].begin(); it2 != face2d_info[intmat_to_string(face2d_nfmat)].end(); ++it2) {
                        std::string key = it2->first;
                        double value = it2->second;
                        // face2d_info[face2d_id][key] = value;
                        // face2d_info[intmat_to_string(face2d_nfmat)][key] = value;
                        if (it == face2d_fulldmat_to_face2d_nfmat.end()) {
                            face2d_fields[key].push_back(value);
                        }
                    }
                    // face2d_fulldmat_to_face2d_nfmat[intmat_to_string(face2d_fulldmat)] = face2d_id;
                    face2d_fulldmat_to_face2d_nfmat[intmat_to_string(face2d_fulldmat)] = intmat_to_string(face2d_nfmat);

                    face2d_fulldmat_to_cone4d_nfmats[intmat_to_string(face2d_fulldmat)].push_back(intmat_to_string(cone4d_nfmat));
                }

                // std::sort(this_face2d_ids.begin(), this_face2d_ids.end());

                // ++cone4d_nfmat_to_mult[face3d_id][vector_to_string(this_face2d_ids)];

                ++cone4d_nfmat_to_mult[intmat_to_string(cone4d_nfmat)];

                // std::string key = std::to_string(face3d_id);
                // if (cone4d_nfmat_to_mult.isMember(key)) {
                //     cone4d_nfmat_to_mult[key].append(inner);
                // } else {
                //     Json::Value outer(Json::arrayValue);
                //     outer.append(inner);
                //     cone4d_nfmat_to_mult[key] = outer;
                // }
            }

            // std::map<int, std::map<std::string, std::vector<double> > > face2d_face3d_info;
            std::map<std::string, std::map<std::string, std::vector<double> > > face2d_face3d_info;
            // std::map<int, std::map<std::string, int> > face2d_nfmat_to_cone4d_nfmats_mult;
            std::map<std::string, std::map<std::string, int> > face2d_nfmat_to_cone4d_nfmats_mult;
            // for (std::map<std::string, int>::iterator it1 = face2d_fulldmat_to_face2d_nfmat.begin(); it1 != face2d_fulldmat_to_face2d_nfmat.end(); it1++) {
            for (std::map<std::string, std::string>::iterator it1 = face2d_fulldmat_to_face2d_nfmat.begin(); it1 != face2d_fulldmat_to_face2d_nfmat.end(); it1++) {
                // int face2d_id = it1->second;
                std::string str_face2d_nfmat = it1->second;
                std::vector<std::string> str_cone4d_nfmats = face2d_fulldmat_to_cone4d_nfmats[it1->first];

                std::map<std::string, double> this_face2d_face3d_info;
                std::vector<std::string> this_cone4d_nfmats;
                for (std::string str_cone4d_nfmat : str_cone4d_nfmats) {
                    std::map<std::string, double> this_face3d_info = face3d_info[str_cone4d_nfmat];
                    for (std::map<std::string, double>::iterator it2 = this_face3d_info.begin(); it2 != this_face3d_info.end(); ++it2) {
                        std::string key = it2->first;
                        double value = it2->second;
                        this_face2d_face3d_info[key] += value;
                    }
                    this_cone4d_nfmats.push_back(str_cone4d_nfmat);
                }
                this_face2d_face3d_info["3FACEN3VERTS"] -= face2d_info[str_face2d_nfmat]["2FACEN2VERTS"];
                this_face2d_face3d_info["3FACEN3POINTS"] -= face2d_info[str_face2d_nfmat]["2FACEN2POINTS"];
                this_face2d_face3d_info["3FACEN2POINTS"] -= face2d_info[str_face2d_nfmat]["2FACEN2POINTS"];
                this_face2d_face3d_info["3FACEN1POINTS"] -= face2d_info[str_face2d_nfmat]["2FACEN1POINTS"];
                this_face2d_face3d_info["3FACEN2FACES"] -= 1;
                this_face2d_face3d_info["3FACEN1FACES"] -= face2d_info[str_face2d_nfmat]["2FACEN1FACES"];

                face2d_face3d_info[str_face2d_nfmat]["COMPLEXN3VERTS"].push_back(this_face2d_face3d_info["3FACEN3VERTS"]);
                face2d_face3d_info[str_face2d_nfmat]["COMPLEXN3POINTS"].push_back(this_face2d_face3d_info["3FACEN3POINTS"]);
                face2d_face3d_info[str_face2d_nfmat]["COMPLEXN2POINTS"].push_back(this_face2d_face3d_info["3FACEN2POINTS"]);
                face2d_face3d_info[str_face2d_nfmat]["COMPLEXN1POINTS"].push_back(this_face2d_face3d_info["3FACEN1POINTS"]);
                face2d_face3d_info[str_face2d_nfmat]["COMPLEX3VOLUME"].push_back(this_face2d_face3d_info["3FACE3VOLUME"]);
                face2d_face3d_info[str_face2d_nfmat]["COMPLEXN2FACES"].push_back(this_face2d_face3d_info["3FACEN2FACES"]);
                face2d_face3d_info[str_face2d_nfmat]["COMPLEXN1FACES"].push_back(this_face2d_face3d_info["3FACEN1FACES"]);

                for (std::map<std::string, std::vector<double> >::iterator it2 = face2d_face3d_info[str_face2d_nfmat].begin(); it2 != face2d_face3d_info[str_face2d_nfmat].end(); ++it2) {
                    std::string key = it2->first;
                    std::vector<double> values = it2->second;
                    for (double value : values) {
                        face2d_fields[key].push_back(value);
                    }
                }

                std::sort(this_cone4d_nfmats.begin(), this_cone4d_nfmats.end());

                ++face2d_nfmat_to_cone4d_nfmats_mult[str_face2d_nfmat][join(this_cone4d_nfmats, ";")];

                // Json::Value inner(Json::arrayValue);
                // for (int face3d_id : face3d_ids) {
                //     inner.append(Json::Value(face3d_id));
                // }
                // if (face2d_nfmat_to_cone4d_nfmats_mult.isMember(key)) {
                //     face2d_nfmat_to_cone4d_nfmats_mult[key].append(inner);
                // } else {
                //     Json::Value outer(Json::arrayValue);
                //     outer.append(inner);
                //     face2d_nfmat_to_cone4d_nfmats_mult[key] = outer;
                // }
            }

            Json::Value json_cone4d_nfmat_to_mult;
            for (std::map<std::string, int>::iterator it1 = cone4d_nfmat_to_mult.begin(); it1 != cone4d_nfmat_to_mult.end(); ++it1) {
                std::string str_cone4d_nfmat = it1->first;
                int mult = it1->second;
                std::string key = str_cone4d_nfmat;
                json_cone4d_nfmat_to_mult[key] = Json::Value(mult);
            }

            Json::Value json_face2d_nfmat_to_cone4d_nfmats_mult;
            // for (std::map<int, std::map<std::string, int> >::iterator it1 = face2d_nfmat_to_cone4d_nfmats_mult.begin(); it1 != face2d_nfmat_to_cone4d_nfmats_mult.end(); ++it1) {
            for (std::map<std::string, std::map<std::string, int> >::iterator it1 = face2d_nfmat_to_cone4d_nfmats_mult.begin(); it1 != face2d_nfmat_to_cone4d_nfmats_mult.end(); ++it1) {
                // int face2d_id = it1->first;
                std::string str_face2d_nfmat = it1->first;
                std::map<std::string, int> cone4d_nfmats_mult = it1->second;
                Json::Value face2d_nfmat_arr(Json::arrayValue);
                for (std::map<std::string, int>::iterator it2 = cone4d_nfmats_mult.begin(); it2 != cone4d_nfmats_mult.end(); ++it2) {
                    std::vector<std::string> cone4d_nfmats = split(it2->first, ";");
                    int mult = it2->second;
                    Json::Value json_cone4d_nfmats_mult(Json::arrayValue);
                    for (std::string cone4d_nfmat : cone4d_nfmats) {
                        json_cone4d_nfmats_mult.append(Json::Value(cone4d_nfmat));
                    }
                    json_cone4d_nfmats_mult.append(Json::Value(mult));
                    face2d_nfmat_arr.append(json_cone4d_nfmats_mult);
                }
                // std::string key = std::to_string(face2d_id);
                std::string key = str_face2d_nfmat;
                json_face2d_nfmat_to_cone4d_nfmats_mult[key] = face2d_nfmat_arr;
            }

            Json::Value index_doc;
            index_doc["POLYID"] = Json::Value(poly_id);

            Json::Value out_doc;
            out_doc["H11"] = Json::Value(h11);
            out_doc["POLYID"] = Json::Value(poly_id);
            out_doc["NVERTS"] = Json::Value(nverts);
            out_doc["4POLYN4VERTS"] = dual_poly.nvertices();
            out_doc["4POLYN4POINTS"] = dual_poly.npoints();
            out_doc["4POLYN3POINTS"] = dual_poly.p_boundary_indices(4).size();
            out_doc["4POLYN2POINTS"] = dual_poly.p_boundary_indices(3).size();
            out_doc["4POLYN1POINTS"] = dual_poly.p_boundary_indices(2).size();
            out_doc["4POLYN3FACES"] = dual_poly.faces(3).size();
            out_doc["4POLYN2FACES"] = dual_poly.faces(2).size();
            out_doc["4POLYN1FACES"] = dual_poly.faces(1).size();
            out_doc["LNNFRTSUM"] = Json::Value(ln_nfrt_sum);
            out_doc["LNNFRTPREDICTSUM"] = Json::Value(ln_nfrt_predict_sum);
            for (std::map<std::string, std::vector<double> >::iterator it = face3d_fields.begin(); it != face3d_fields.end(); ++it) {
                std::string key = it->first;
                std::vector<double> values = it->second;
                double min = *std::min_element(values.begin(), values.end());
                double max = *std::max_element(values.begin(), values.end());
                double sum = std::accumulate(values.begin(), values.end(), 0.0);
                double mean = sum / values.size();
                std::vector<double> diffsq_values;
                for (int i = 0; i < values.size(); i++) {
                    diffsq_values.push_back(std::pow(values[i] - mean, 2));
                }
                double diffsq_sum = std::accumulate(diffsq_values.begin(), diffsq_values.end(), 0.0);
                double stdev = std::sqrt(diffsq_sum / (diffsq_values.size() - 1));
                out_doc[key]["MIN"] = min;
                out_doc[key]["MAX"] = max;
                out_doc[key]["MEAN"] = mean;
                out_doc[key]["STDEV"] = stdev;
            }
            for (std::map<std::string, std::vector<double> >::iterator it = face2d_fields.begin(); it != face2d_fields.end(); ++it) {
                std::string key = it->first;
                std::vector<double> values = it->second;
                double min = *std::min_element(values.begin(), values.end());
                double max = *std::max_element(values.begin(), values.end());
                double sum = std::accumulate(values.begin(), values.end(), 0.0);
                double mean = sum / values.size();
                std::vector<double> diffsq_values;
                for (int i = 0; i < values.size(); i++) {
                    diffsq_values.push_back(std::pow(values[i] - mean, 2));
                }
                double diffsq_sum = std::accumulate(diffsq_values.begin(), diffsq_values.end(), 0.0);
                double stdev = std::sqrt(diffsq_sum / (diffsq_values.size() - 1));
                out_doc[key]["MIN"] = min;
                out_doc[key]["MAX"] = max;
                out_doc[key]["MEAN"] = mean;
                out_doc[key]["STDEV"] = stdev;
            }

            if (poly_doc.isMember("NFSRT")) {
                int poly_triang = poly_doc["NFSRT"].asInt();
                out_doc["NFSRT"] = Json::Value(poly_triang);
            }

            // out_doc["CONETOFACE"] = json_cone4d_nfmat_to_mult;
            // out_doc["FACETOCONE"] = json_face2d_nfmat_to_cone4d_nfmats_mult;

            std::cout << "set OVERCOUNT ";
            writer->write(index_doc, &std::cout);
            std::cout << " ";
            writer->write(out_doc, &std::cout);
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}