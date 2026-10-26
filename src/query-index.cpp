/*
 * This is free software: you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the
 * Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This software is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 * See the GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <iostream>
#include <utility>
#include <ring_fixed_p.hpp>
#include <chrono>
#include <triple_pattern.hpp>
#include <ltj_algorithm.hpp>
#include <utils.hpp>

using namespace std;

//#include<chrono>
//#include<ctime>

using namespace std::chrono;

bool get_file_content(string filename, vector<string> & vector_of_strings)
{
    // Open the File
    ifstream in(filename.c_str());
    // Check if object is valid
    if(!in)
    {
        cerr << "Cannot open the File : " << filename << endl;
        return false;
    }
    string str;
    // Read the next line from File until it reaches the end.
    while (getline(in, str))
    {
        // Line contains string of length > 0 then save it in vector
        if(str.size() > 0)
            vector_of_strings.push_back(str);
    }
    //Close The File
    in.close();
    return true;
}

std::string ltrim(const std::string &s)
{
    size_t start = s.find_first_not_of(' ');
    return (start == std::string::npos) ? "" : s.substr(start);
}

std::string rtrim(const std::string &s)
{
    size_t end = s.find_last_not_of(' ');
    return (end == std::string::npos) ? "" : s.substr(0, end + 1);
}

std::string trim(const std::string &s) {
    return rtrim(ltrim(s));
}

std::vector<std::string> tokenizer(const std::string &input, const char &delimiter){
    std::stringstream stream(input);
    std::string token;
    std::vector<std::string> res;
    while(getline(stream, token, delimiter)){
        res.emplace_back(trim(token));
    }
    return res;
}

bool is_variable(string & s)
{
    return (s.at(0) == '?');
}

uint8_t get_variable(string &s, std::unordered_map<std::string, uint8_t> &hash_table_vars){
    auto var = s.substr(1);
    auto it = hash_table_vars.find(var);
    if(it == hash_table_vars.end()){
        uint8_t id = hash_table_vars.size();
        hash_table_vars.insert({var, id });
        return id;
    }else{
        return it->second;
    }
}

uint64_t get_constant(string &s){
    return std::stoull(s);
}

ring::triple_pattern get_triple(string & s, std::unordered_map<std::string, uint8_t> &hash_table_vars) {
    vector<string> terms = tokenizer(s, ' ');

    ring::triple_pattern triple;
    if(is_variable(terms[0])){
        triple.var_s(get_variable(terms[0], hash_table_vars));
    }else{
        triple.const_s(get_constant(terms[0]));
    }
    if(is_variable(terms[1])){
        triple.var_p(get_variable(terms[1], hash_table_vars));
    }else{
        triple.const_p(get_constant(terms[1]));
    }
    if(is_variable(terms[2])){
        triple.var_o(get_variable(terms[2], hash_table_vars));
    }else{
        triple.const_o(get_constant(terms[2]));
    }
    return triple;
}


typedef ring::ring_fixed_p ring_type;
typedef ring::ltj_iterator<ring_type, uint8_t, uint64_t> iter_type;

template<class veo_type>
void run_queries(ring_type &graph, vector<string> &queries, const uint64_t limit, const uint64_t timeout) {

    typedef ring::ltj_algorithm<ring_type, uint8_t, uint64_t, veo_type> algorithm_type;
    typedef std::vector<typename algorithm_type::tuple_type> results_type;

    high_resolution_clock::time_point start, stop;
    uint64_t nQ = 0;
    for (string& query_string : queries) {

        std::unordered_map<std::string, uint8_t> hash_table_vars;
        std::vector<ring::triple_pattern> query;
        vector<string> tokens_query = tokenizer(query_string, '.');
        for (string& token : tokens_query) {
            auto triple_pattern = get_triple(token, hash_table_vars);
            query.push_back(triple_pattern);
        }

        ring::veo::stats::reset();
        start = high_resolution_clock::now();
        algorithm_type ltj(&query, &graph);
        results_type res;
        ltj.join_opt(res, limit, timeout);
        stop = high_resolution_clock::now();
        auto ns = duration_cast<nanoseconds>(stop - start).count();

        cout << nQ << ";" << res.size() << ";" << ns;
#ifdef VEO_STATS
        cout << ";" << ring::veo::stats::calls() << ";" << ring::veo::stats::ns();
#endif
        cout << endl;
        nQ++;
    }
}

template<uint64_t k>
bool run_k(const string &veo, ring_type &graph, vector<string> &queries, const uint64_t limit, const uint64_t timeout) {
    using namespace ring::veo;
    const string K = std::to_string(k);
    if (veo == "global-int" + K) run_queries<veo_global<iter_type, estimator_intersect<k, false>>>(graph, queries, limit, timeout);
    else if (veo == "global-guided" + K) run_queries<veo_global<iter_type, estimator_intersect<k, true>>>(graph, queries, limit, timeout);
    else if (veo == "adaptive-int" + K) run_queries<veo_adaptive<iter_type, estimator_intersect<k, false>>>(graph, queries, limit, timeout);
    else if (veo == "adaptive-guided" + K) run_queries<veo_adaptive<iter_type, estimator_intersect<k, true>>>(graph, queries, limit, timeout);
    else return false;
    return true;
}

int main(int argc, char* argv[])
{
    if (argc < 3) {
        cerr << "Usage: " << argv[0] << " <index> <queries> [limit=1000 (0: no limit)] [veo=adaptive-guided3] [timeout=600 s]" << endl;
        cerr << "  veo: global | adaptive | {global,adaptive}-{int,guided}{1,2,3,4}" << endl;
        cerr << "    global:   order fixed before running LTJ (original RingFixedP order with global);" << endl;
        cerr << "    adaptive: next variable chosen after each binding" << endl;
        cerr << "    int<k>:   refined estimator, 2^k equal parts of [0,U];  guided<k>: 2^k quantiles of the smallest interval" << endl;
        return 1;
    }
    uint64_t limit = (argc > 3) ? std::stoull(argv[3]) : 1000;
    string veo = (argc > 4) ? argv[4] : "adaptive-guided3";
    uint64_t timeout = (argc > 5) ? std::stoull(argv[5]) : 600;

    vector<string> queries;
    if (!get_file_content(argv[2], queries)) return 1;

    ring_type graph;
    cout << " Loading the index..."; fflush(stdout);
    sdsl::load_from_file(graph, string(argv[1]));
    cout << endl << " Index loaded " << sdsl::size_in_bytes(graph) << " bytes" << endl;

    using namespace ring::veo;
    if (veo == "global") run_queries<veo_global<iter_type, estimator_size>>(graph, queries, limit, timeout);
    else if (veo == "adaptive") run_queries<veo_adaptive<iter_type, estimator_size>>(graph, queries, limit, timeout);
    else if (!run_k<1>(veo, graph, queries, limit, timeout) && !run_k<2>(veo, graph, queries, limit, timeout)
             && !run_k<3>(veo, graph, queries, limit, timeout) && !run_k<4>(veo, graph, queries, limit, timeout)
) {
        cerr << "Unknown veo: " << veo << endl;
        return 1;
    }
    return 0;
}
