#include <iostream>
#include <chrono>
#include <ctime>
#include <vector>
#include <sstream>
#include <fstream>

// Security define
#define AES_SECURITY 80
#define MR_PAIRING_SSP
#include "pairing_1.h"
//#include "ecn.h"

using namespace std;
//extern "C" {
//#include "miracl.h"
//#include "mirdef.h"
//}

struct KeyParams
{
    Big q;
    G1 P;
};

struct Pi_list
{
    GT M;
    GT N;
    GT R1;
    vector<G1> U_array;
    Big e;
};

struct Pi2_list
{
    GT M;
    GT N;
    Big e;
};

struct Pi3_list
{
    GT Mi;
    GT Ni;
    Big ei;
};

struct Pi3
{
    GT pi;
    Pi3_list pi3_list;
};

struct Sigma
{
    GT p;
    Big r0;
    Pi_list pi_list;
};

pair<Big, G1> KeyGen(KeyParams& params, PFC& pfc) {
    Big x;
    pfc.random(x);
    G1 y = pfc.mult(params.P, x);
    return { x, y };
}


void add_to_hash(sha256& sh, const G1& x)
{
    Big a, X, Y;
    int i, m;
    x.g.get(X, Y);
    a = X;
    while (a > 0)
    {
        m = a % 256;
        shs256_process(&sh, m);
        a /= 256;
    }
    a = Y;
    while (a > 0)
    {
        m = a % 256;
        shs256_process(&sh, m);
        a /= 256;
    }
}

void add_to_hash(sha256& sh, const Big& x)
{
    int m;
    Big a = x;
    while (a > 0)
    {
        m = a % 256;
        shs256_process(&sh, m);
        a /= 256;
    }
}

G1 H0(KeyParams& params, PFC& pfc, const Big& bit, const Big& r0, char* m, const vector<G1>& R) {
    // Convert input parameters to a single string
    G1 G1_result;
    Big temp_m = pfc.hash_to_group(m);

    sha256 sh;
    shs256_init(&sh);
    add_to_hash(sh, bit);
    add_to_hash(sh, r0);
    add_to_hash(sh, temp_m);
    for (const G1& r : R) {
        add_to_hash(sh, r);
    }

    char hash_result[32];
    shs256_hash(&sh, hash_result);

    pfc.hash_and_map(G1_result, hash_result);

    return G1_result;
}


//G1 H0(KeyParams& params, PFC& pfc, Big& bit, Big& r0, char* m, const vector<G1>& pk_list) {
//    // Convert Big parameters to strings
//    G1 result;
//    stringstream ss;
//    //Big temp_m = pfc.hash_to_group(m);
//    ss << bit << r0 << m;
//    for (const auto& r : pk_list) {
//        Big x, y;
//        r.g.getxy(x, y);
//        ss << x << y;
//    }
//
//    // Convert combined string to a char array
//    
//    string str = ss.str();
//
//    
//
//    //result = HashToG1(pfc, str);
//    const char* str_char = str.c_str();
//    int str_len = strlen(str_char);
//    int len = str.length();
//
//    //char* combined_char = malloc(strlen(str_char) + 1)
//
//    char* combined_char = new char[len + 1];
//    //vector<string> char_array(combined_char);
//    //int com_len = sizeof(combined_char);
//    strncpy_s(combined_char, len + 1, str_char, len);
//    //strncpy_s(char_array, str.c_str(), len);
//    //strncpy_s(combined_char, ss.str().c_str(), 1001);
//    //char* char_array = new char[str.length() + 1];
//    ////std::strcpy(char_array, str.c_str());
//    //strcpy_s(combined_char, len + 1, str.c_str());
//    //strcat_s(combined_char, len+1, str_char);
//    //combined_char[len-1] = '\0'; // Ensure null-termination
//    // 
//    //int com_len = strlen(combined_char);
//    // Hash and map to G1
//    
//    pfc.hash_and_map(result, combined_char);
//    //_ASSERTE(_CrtCheckMemory());
//    delete[] combined_char;
//    delete[] str_char;
//    return result;
//}

Big H(KeyParams& params, PFC& pfc, char* m, const GT& M, const GT& N, const GT& R1, const GT& p, const G1& Ui) {
    Big temp_m = pfc.hash_to_group(m);
    pfc.start_hash();
    pfc.add_to_hash(temp_m);
    pfc.add_to_hash(M);
    pfc.add_to_hash(N);
    pfc.add_to_hash(R1);
    pfc.add_to_hash(p);
    pfc.add_to_hash(Ui);
    Big result = pfc.finish_hash_to_group();
    return result;
}


Sigma sign(KeyParams& params, PFC& pfc, vector<G1>& pk_list, Big& sk, int k, char* m) {
    Big r0 = rand(params.q);
    G1 P = params.P;
    int len = pk_list.size();
    //Big r0;
    //pfc.random(r0);
        
    Big bit0 = 0;
    Big bit1 = 0;
    G1 u0 = H0(params, pfc, bit0, r0, m, pk_list);
    //G1 u00 = H0(params, pfc, bit0, r0, m, pk_list);
    G1 u1= H0(params, pfc, bit1,r0, m, pk_list);

    /*if (u0 == u00)
    {
        cout << "this is eq!" << endl;
    }*/

    GT p = pfc.pairing(u1, u0);
    p = pfc.power(p, sk);

    Big d, r1;
    pfc.random(d);
    pfc.random(r1);
    GT M, N, R1;
    M = pfc.pairing(P, P);
    M = pfc.power(M, d);
    N = pfc.pairing(u1, u0);
    N = pfc.power(N, d);
    //R1 = pfc.pairing(p, P);
    R1 = pfc.power(p, r1);

    vector<G1> U_array(len);
    vector<Big> h_array(len);
    G1 sum_R;
    sum_R.g.clear();


    for (int i = 0; i < len; i++)
    {
        if (i == k)
        {
            continue;
        }
        pfc.random(U_array[i]);
        h_array[i] = H(params, pfc, m, M, N, R1, p, U_array[i]);
        G1 sum_temp = U_array[i] + pfc.mult(pk_list[i], h_array[i]);
        G1 pkk_temp = -pfc.mult(pk_list[k], h_array[i]);
        
        sum_R = sum_R + sum_temp + pkk_temp;

            
    }
    G1 sum_R_invers = -sum_R;
    U_array[k] = pfc.mult(pk_list[k], r1) + sum_R_invers;
    
    h_array[k] = H(params, pfc, m, M, N, R1, p, U_array[k]);
    Big sum_h = 0;
    for (int i = 0; i < len; i++)
    {
        sum_h += h_array[i];
    }
    Big e = d - (sk * (sum_h + r1));
    Pi_list pi_list;
    pi_list.M = M;
    pi_list.N = N;
    pi_list.R1 = R1;
    pi_list.U_array = U_array;
    pi_list.e = e;
    return { p, r0, pi_list };
    
}

bool verify(KeyParams& params, PFC& pfc, vector<G1>& pk_list, char* m, Sigma& sigma) {
    G1 P = params.P;
    GT p = sigma.p;
    Big r0 = sigma.r0;
    Pi_list pi_list = sigma.pi_list;
    GT M = pi_list.M;
    GT N = pi_list.N;
    GT R1 = pi_list.R1;
    vector<G1> U_array = pi_list.U_array;
    Big e = pi_list.e;

    G1 sum_R;
    sum_R.g.clear();
    Big sum_h = 0;
    for (int i = 0; i < pk_list.size(); i++)
    {
        Big hi = H(params, pfc, m, M, N, R1, p, U_array[i]);
        sum_h += hi;
        sum_R = sum_R + U_array[i] + pfc.mult(pk_list[i], hi);
    }

    Big bit0 = 0;
    Big bit1 = 0;
    G1 u0 = H0(params, pfc, bit0, r0, m, pk_list);
    G1 u1 = H0(params, pfc, bit1, r0, m, pk_list);

    GT righteq11 = pfc.pairing(P, P);
    righteq11 = pfc.power(righteq11, e);
    GT righteq12 = pfc.pairing(P, sum_R);
    GT righteq1 = righteq11 * righteq12;
    if (M != righteq1)
    {
        cout << "equation 1 verify error!!!" << endl;
        return false;
    }

    GT righteq21 = pfc.pairing(u1, u0);
    righteq21 = pfc.power(righteq21, e);
    GT righteq22 = pfc.power(p, sum_h);
    GT righteq2 = righteq21 * R1 * righteq22;
    if (N != righteq2)
    {
        cout << "equation 2 verify error!!!" << endl;
        return false;
    }
    return true;
}

Pi2_list Vk_prove(KeyParams& params, PFC& pfc, vector<G1>& pk_list, char* m, Sigma& sigma, Big sk) {
    G1 P = params.P;
    GT p = sigma.p;
    Big r0 = sigma.r0;

    Big bit0 = 0;
    Big bit1 = 0;
    G1 u0 = H0(params, pfc, bit0, r0, m, pk_list);
    G1 u1 = H0(params, pfc, bit1, r0, m, pk_list);

    Big d;
    pfc.random(d);
    GT M, N;
    M = pfc.pairing(P, P);
    M = pfc.power(M, d);
    N = pfc.pairing(u1, u0);
    N = pfc.power(N, d);
    pfc.start_hash();
    pfc.add_to_hash(M);
    pfc.add_to_hash(N);
    pfc.add_to_hash(p);
    Big hk = pfc.finish_hash_to_group();

    Big e = d - hk * sk;

    Pi2_list pi2_list;
    pi2_list.e = e;
    pi2_list.M = M;
    pi2_list.N = N;

    return pi2_list;

}

// TRC verify
bool TRC_trace1(KeyParams& params, PFC& pfc, vector<G1>& pk_list, char* m, Sigma& sigma, Pi2_list& pi2_list, int k) {
    G1 P = params.P;
    GT p = sigma.p;
    Big r0 = sigma.r0;
    Big e = pi2_list.e;
    GT M = pi2_list.M;
    GT N = pi2_list.N;
    pfc.start_hash();
    pfc.add_to_hash(M);
    pfc.add_to_hash(N);
    pfc.add_to_hash(p);
    Big hk = pfc.finish_hash_to_group();

    Big bit0 = 0;
    Big bit1 = 0;
    G1 u0 = H0(params, pfc, bit0, r0, m, pk_list);
    G1 u1 = H0(params, pfc, bit1, r0, m, pk_list);

    int flag1 = 0, flag2 = 0;

    // verify pi2
    GT righteq11 = pfc.pairing(P, P);
    righteq11 = pfc.power(righteq11, e);
    GT righteq12 = pfc.pairing(P, pk_list[k]);
    righteq12 = pfc.power(righteq12, hk);
    GT righteq1 = righteq11 * righteq12;

    if (M == righteq1)
    {
        flag1 = 1;
    }

    GT righteq21 = pfc.pairing(u1, u0);
    righteq21 = pfc.power(righteq21, e);
    GT righteq22 = pfc.power(p, hk);
    GT righteq2 = righteq21 * righteq22;
    if (N == righteq2)
    {
        flag2 = 1;
    }

    if (flag1 == 1 && flag2 == 1)
    {
        return true;
    }
    return false;
}

Pi3 disavowal(KeyParams& params, PFC& pfc, vector<G1>& pk_list, char* m, Sigma& sigma, Big ski) {
    G1 P = params.P;

    Big r0 = sigma.r0;

    Big bit0 = 0;
    Big bit1 = 0;
    G1 u0 = H0(params, pfc, bit0, r0, m, pk_list);
    G1 u1 = H0(params, pfc, bit1, r0, m, pk_list);

    GT pi = pfc.pairing(u1, u0);
    pi = pfc.power(pi, ski);

    Big di;
    pfc.random(di);
    GT Mi, Ni;
    Mi = pfc.pairing(P, P);
    Mi = pfc.power(Mi, di);
    Ni = pfc.pairing(u1, u0);
    Ni = pfc.power(Ni, di);
    pfc.start_hash();
    pfc.add_to_hash(Mi);
    pfc.add_to_hash(Ni);
    pfc.add_to_hash(pi);
    Big hki = pfc.finish_hash_to_group();
    Big ei = di - hki * ski;

    Pi3_list pi3_list;
    pi3_list.ei = ei;
    pi3_list.Mi = Mi;
    pi3_list.Ni = Ni;

    return { pi, pi3_list };
}

bool TRC_trace2(KeyParams& params, PFC& pfc, vector<G1>& pk_list, char* m, Sigma& sigma, Pi3& pi3, int i) {
    G1 P = params.P;
    GT p = sigma.p;
    Big r0 = sigma.r0;
    Pi3_list pi3_list = pi3.pi3_list;
    GT pi = pi3.pi;
    Big ei = pi3_list.ei;
    GT Mi = pi3_list.Mi;
    GT Ni = pi3_list.Ni;
    pfc.start_hash();
    pfc.add_to_hash(Mi);
    pfc.add_to_hash(Ni);
    pfc.add_to_hash(pi);
    Big hki = pfc.finish_hash_to_group();

    Big bit0 = 0;
    Big bit1 = 0;
    G1 u0 = H0(params, pfc, bit0, r0, m, pk_list);
    G1 u1 = H0(params, pfc, bit1, r0, m, pk_list);

    int flag1 = 0, flag2 = 0;

    // verify pi2
    GT righteq11 = pfc.pairing(P, P);
    righteq11 = pfc.power(righteq11, ei);
    GT righteq12 = pfc.pairing(P, pk_list[i]);
    righteq12 = pfc.power(righteq12, hki);
    GT righteq1 = righteq11 * righteq12;

    if (Mi == righteq1)
    {
        flag1 = 1;
    }

    GT righteq21 = pfc.pairing(u1, u0);
    righteq21 = pfc.power(righteq21, ei);
    GT righteq22 = pfc.power(pi, hki);
    GT righteq2 = righteq21 * righteq22;
    if (Ni == righteq2)
    {
        flag2 = 1;
    }
    if (pi != p)
    {
        if (flag1 == 1 && flag2 == 1)
        {
            return true;
        }
    }
    
    return false;
}

void test_trial(KeyParams& params, PFC& pfc) {
    vector<double> sign_time_list;
    vector<double> trace1_time_list;
    vector<double> trace2_time_list;
    vector<double> verify_time_list;
    vector<double> entire_time_list;
    vector<double> confirmation_time_list;
    vector<double> disavowal_time_list;

    cout << "power\tsign\tverify\tconfirmation\ttrace1\tdisavowal\ttrace2\tentire(ms)" << endl;

    for (int power = 0; power < 6; power++) {
        double sign_time_sum = 0;
        double trace1_time_sum = 0;
        double trace2_time_sum = 0;
        double verify_time_sum = 0;
        double entire_time_sum = 0;
        double confirmation_time_sum = 0;
        double disavowal_time_sum = 0;
        int power_of_2 = power + 1;
        int PK_num = (int)pow(2, power_of_2);
        int time_trail = 10;

        vector<G1> fake_Pid(PK_num);
        vector<Big> fake_Ssk(PK_num);

        pair<Big, G1> psk;
        for (int j = 0; j < PK_num; j++) {
            psk = KeyGen(params, pfc);
            fake_Pid[j] = psk.second;
            fake_Ssk[j] = psk.first;
        }

        time_t seed;
        time(&seed);
        irand((long)seed);


        for (int ii = 0; ii < time_trail; ii++) {
            clock_t start_time = clock();
            int random_position = rand(params.q) % PK_num;
            
            psk = KeyGen(params, pfc);
            fake_Pid[random_position] = psk.second;
            Big ssk = psk.first;

           
            clock_t sign_start_time = clock();
            
            Sigma sigma = sign(params, pfc, fake_Pid, ssk, random_position, (char*)"I am a girl");
            clock_t sign_end_time = clock();
            sign_time_sum += ((double)sign_end_time - (double)sign_start_time) / CLOCKS_PER_SEC;

            /*clock_t confirmation_start_time = clock();

            Pi2_list pi2 = Vk_prove(params, pfc, fake_Pid, (char*)"I am a girl", sigma, ssk);
            clock_t confirmation_end_time = clock();
            confirmation_time_sum += ((double)confirmation_end_time - (double)confirmation_start_time) / CLOCKS_PER_SEC;*/

            //vector<Pi3> pi3_array(PK_num);
            //clock_t disavowal_start_time = clock();

            //for (int iii = 0; iii < PK_num; iii++)
            //{
            //    pi3_array[iii] = disavowal(params, pfc, fake_Pid, (char*)"I am a girl", sigma, fake_Ssk[iii]);
            //}
            //
            //clock_t disavowal_end_time = clock();
            //disavowal_time_sum += ((double)disavowal_end_time - (double)disavowal_start_time) / CLOCKS_PER_SEC;

            ///*clock_t trace1_start_time = clock();
            //if (!TRC_trace1(params, pfc, fake_Pid, (char*)"I am a girl", sigma, pi2, random_position)) {
            //    cout << "Confirmation, TRC trace failed" << endl;
            //}
            //
            //clock_t trace1_end_time = clock();
            //trace1_time_sum += ((double)trace1_end_time - (double)trace1_start_time) / CLOCKS_PER_SEC;*/

            //clock_t trace2_start_time = clock();
            //for (int iv = 0; iv < PK_num; iv++)
            //{
            //    if (!TRC_trace2(params, pfc, fake_Pid, (char*)"I am a girl", sigma, pi3_array[iv], iv)) {
            //        //cout << "TRC trace success! The malicious user index is "<< iv << endl;
            //    }
            //}
            //
            //clock_t trace2_end_time = clock();
            //trace2_time_sum += ((double)trace2_end_time - (double)trace2_start_time) / CLOCKS_PER_SEC;

            clock_t verify_start_time = clock();
            if (!verify(params, pfc, fake_Pid, (char*)"I am a girl", sigma)) {
                cout << "signature verify failed" << endl;
            }
            clock_t verify_end_time = clock();
            verify_time_sum += ((double)verify_end_time - (double)verify_start_time) / CLOCKS_PER_SEC;

            clock_t end_time = clock();
            entire_time_sum += ((double)end_time - (double)start_time) / CLOCKS_PER_SEC;
        }

        cout << power + 1 << "\t" << (sign_time_sum / time_trail) * 1000 << "\t" << (verify_time_sum / time_trail) * 1000
            << "\t" << (confirmation_time_sum / time_trail) * 1000 << "\t" << (trace1_time_sum / time_trail) * 1000  
            << "\t" << (disavowal_time_sum / time_trail) * 1000 << "\t" << (trace2_time_sum / time_trail) * 1000 
            << "\t" << (entire_time_sum / time_trail) * 1000 << endl;

        sign_time_list.push_back((sign_time_sum / time_trail) * 1000);
        verify_time_list.push_back((verify_time_sum / time_trail) * 1000);
        confirmation_time_list.push_back((confirmation_time_sum / time_trail) * 1000);
        trace1_time_list.push_back((trace1_time_sum / time_trail) * 1000);
        disavowal_time_list.push_back((disavowal_time_sum / time_trail) * 1000);
        trace2_time_list.push_back((trace2_time_sum / time_trail) * 1000);
        entire_time_list.push_back((entire_time_sum / time_trail) * 1000);
    }

    // Get current date
    std::time_t t = std::time(nullptr);
    std::tm now;
    localtime_s(&now, &t);

    // Format date as YYYY-MM-DD
    std::ostringstream date_stream;
    date_stream << (now.tm_year + 1900) << '-'
        << (now.tm_mon + 1) << '-'
        << now.tm_mday << '-' << now.tm_hour << '-' << now.tm_min << '-';
    std::string date_str = date_stream.str();

    // Create filename with date
    std::string filename = date_str + "CARS Ring Signature Time Analysis .txt";

    // 指定文件路径
    std::string directory = "D:\\VS2019\\repo\\2024-6-10-CARS\\test_time\\";
    //std::string filename = "output.txt";
    std::string filepath = directory + filename;

    // Writing to file
    std::ofstream text_file;
    text_file.open(filepath);

    // Writing to file
    //std::ofstream text_file(filename);
    if (!text_file.is_open()) {
        std::cerr << "Failed to open file: " << filename << std::endl;
        //return 1;
    }

    // writing to file
    //ofstream text_file("DualRing Ring Signature Time Analysis.txt");
    text_file << "2^n\tSign\tVerify\tConfir\tTrace1\tDisa\tTrace2\n";
    for (size_t i = 0; i < sign_time_list.size(); i++) {
        text_file << i + 1 << "\t" << sign_time_list[i] << "\t" << verify_time_list[i] << "\t" << confirmation_time_list[i]
            << "\t" << trace1_time_list[i] << "\t" << disavowal_time_list[i] << "\t" << trace2_time_list[i] << "\t" << endl;
    }
    text_file.close();
}


int main() {
    // Security level for 80 bits
    PFC pfc(80);
    int trials = 1000;
    Big a, b;
    Big q, x_str;
    G1 P, Q, P1;
    pfc.random(P);

    q = pfc.order();
    // TRA secret and public key
    pfc.random(x_str);
    G1 y_str = pfc.mult(P, x_str);

    KeyParams params;
    params.q = q;
    params.P = P;

    test_trial(params, pfc);

    return 0;
}
