#include <bits/stdc++.h>
#include <sys/time.h>
#include <atcoder/dsu>

using namespace std;

const vector<tuple<int,int,int>> d3 = {{0, 0, 1}, {0, 1, 0}, {1, 0, 0}, {0, 0, -1}, {0, -1, 0}, {-1, 0, 0}}; 
vector<vector<vector<vector<int>>>> vaild_pos;
vector<vector<tuple<int,int,int>>> outercells(2);
vector<vector<tuple<int,int,int>>> out8(2, vector<tuple<int,int,int>>(8));

struct rand_generator {
    random_device seed_gen;
    mt19937 engine;
    mt19937_64 engine64;
    rand_generator() : engine(seed_gen()), engine64(seed_gen()) {}
    int rand(int mod) {
        return engine() % mod;
    }
    long long randll(long long mod) {
        return engine64() % mod;
    }
} Ryuka;

struct timer {
    double global_start;
    double gettime() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec * 1e-6;
    }
    void init() {
        global_start = gettime();
    }
    double elapsed() {
        return gettime() - global_start;
    }
} Toki;

struct logger {
    bool enable;
    logger() {
#if defined(ONLINE_JUDGE)
        enable = false;
#else
        enable = true;
#endif
    }
    void info(string str) {
        if(enable) {
            cerr << "[INFO] \t" << str << "\n";
        }
    }
    void raw(string str) {
        if(enable) {
            cerr << str << "\n";
        }
    }
} Sera;

struct Input {
    int D;
    vector<vector<vector<int>>> f, r;
    vector<vector<string>> _f, _r;
    
    void read() {
        cin >> D;
        _f.resize(2, vector<string>(D));
        _r.resize(2, vector<string>(D));
        f.resize(2, vector(D, vector<int>(D)));
        r.resize(2, vector(D, vector<int>(D)));
        for(int i = 0; i < 2; i++) {
            for(int j = 0; j < D; j++) {
                cin >> _f[i][j];
            }
            for(int j = 0; j < D; j++) {
                cin >> _r[i][j];
            }
        }
        for(int i = 0; i < 2; i++) {
            for(int j = 0; j < D; j++) {
                for(int k = 0; k < D; k++) {
                    f[i][j][k] = _f[i][j][k] - '0';
                    r[i][j][k] = _r[i][j][k] - '0';
                }
            }
        }
    }
} inputs;
struct hash {
    vector<vector<vector<long long>>> table;
    void init() {
        table.resize(inputs.D, vector(inputs.D, vector<long long>(inputs.D)));
        for(int x = 0; x < inputs.D; x++) {
            for(int y = 0; y < inputs.D; y++) {
                for(int z = 0; z < inputs.D; z++) {
                    table[x][y][z] = Ryuka.randll(1LL<<60);
                }
            }
        }
    }
    long long get_hash(const vector<tuple<int,int,int>>& v) const {
        long long res = 0;
        for(const auto& [x, y, z] : v) res ^= table[x][y][z];
        return res;
    }
} FunaQ;


struct Output {
    int n;
    vector<vector<vector<vector<int>>>> ans;
    Output() { 
        ans.resize(2, vector(inputs.D, vector(inputs.D, vector<int>(inputs.D, 0)))); 
    }
    Output(int n, vector<vector<vector<vector<int>>>> ans) : n(n), ans(ans) { ; }
    void write() {
        cout << n << endl;
        for(int i = 0; i < 2; i++) {
            for(int x = 0; x < inputs.D; x++) {
                for(int y = 0; y < inputs.D; y++) {
                    for(int z = 0; z < inputs.D; z++) {
                        cout << ans[i][x][y][z] << " ";
                    }
                }
            }
            cout << endl;
        }
    }
};

struct Record {
    int base;
    int initial_points;
    long long score;
    vector<vector<tuple<int,int,int>>> initial_pos;
    Output outputs;
    Record() {
        initial_pos.resize(2);
    }
    Record(int base, int initial_points, long long score, vector<vector<tuple<int,int,int>>> initial_pos) 
        : base(base), initial_points(initial_points), score(score), initial_pos(initial_pos) {
    }
};

namespace Utils {

    vector<tuple<int,int,int>> normalize(vector<tuple<int,int,int>> v) {
        int min_x = 1e9, min_y = 1e9, min_z = 1e9;
        for(const auto &[x, y, z] : v) {
            min_x = min(min_x, x);
            min_y = min(min_y, y);
            min_z = min(min_z, z);
        }
        for(auto &[x, y, z] : v) {
            x -= min_x;
            y -= min_y;
            z -= min_z;
        }
        return v;
    }

    set<pair<long long, long long>> differ;
    int hash_counter = 0, hash_total = 0;

    bool is_same(const vector<tuple<int,int,int>>& a, const vector<tuple<int,int,int>>& b, bool use_hash = true) {
        if(a.size() != b.size()) {
            return false;
        }
        if(inputs.D >= 8) use_hash = false;
        auto b1 = normalize(a);
        auto b2 = normalize(b);
        long long hashA = FunaQ.get_hash(b1);
        long long hashB = FunaQ.get_hash(b2);
        hash_total++;
        if(use_hash && differ.count({hashA, hashB})) {
            hash_counter++;
            return false;
        }
        int max_x = 0, max_y = 0, max_z = 0;
        for(auto &[x, y, z]: b2) {
            max_x = max(max_x, x);
            max_y = max(max_y, y);
            max_z = max(max_z, z);
        }
        sort(b1.begin(), b1.end());
        for(int i = 0; i < 6; i++) {
            for(int _ = 0; _ < 4; _++) {
                sort(b2.begin(), b2.end());
                if(b1 == b2) {
                    return true;
                }
                for(auto &[x, y, _] : b2) {
                    int t = x;
                    x = max_y - y;
                    y = t;
                }
                swap(max_x, max_y);
            }
            if((i & 1) != 0) {
                for(auto &[_, y, z] : b2) {
                    int t = y;
                    y = max_z - z;
                    z = t;
                }
                swap(max_y, max_z);
            } else {
                for(auto &[x, _, z] : b2) {
                    int t = z;
                    z = max_x - x;
                    x = t;
                }
                swap(max_x, max_z);
            }
        }
        if(use_hash) {
            differ.insert({hashA, hashB});
            differ.insert({hashB, hashA});
        }
        return false;
    }

    long long compute_score(const Output& outputs) {
        const auto& n = outputs.n;
        const auto& ans = outputs.ans;
        vector pos(2, vector<vector<tuple<int,int,int>>>(n));
        vector visited(vector(2, vector(inputs.D, vector(inputs.D, vector<bool>(inputs.D, false)))));
        for(int i = 0; i < 2; i++) {
            vector f(inputs.D, vector<int>(inputs.D, 0));
            vector r(inputs.D, vector<int>(inputs.D, 0));
            for(int x = 0; x < inputs.D; x++) {
                for(int y = 0; y < inputs.D; y++) {
                    for(int z = 0; z < inputs.D; z++) {
                        int id = ans[i][x][y][z];
                        if(id != 0) {
                            f[z][x] = 1;
                            r[z][y] = 1;
                            pos[i][id-1].push_back({x, y, z});
                            if(pos[i][id-1].size() == 1) {
                                visited[i][x][y][z] = true;
                                stack<tuple<int,int,int>> stk;
                                stk.push({x, y, z});
                                while(!stk.empty()) {
                                    auto [x, y, z] = stk.top();
                                    stk.pop();
                                    for(const auto &[dx, dy, dz] : d3) {
                                        int x2 = x + dx;
                                        int y2 = y + dy;
                                        int z2 = z + dz;
                                        if (x2 < inputs.D && x2 >= 0 &&
                                            y2 < inputs.D && y2 >= 0 &&
                                            z2 < inputs.D && z2 >= 0 &&
                                            ans[i][x2][y2][z2] == id &&
                                            !visited[i][x2][y2][z2]
                                        )
                                        {
                                            visited[i][x2][y2][z2] = true;
                                            stk.push({x2, y2, z2});
                                        }
                                    }
                                }
                            } else if(!visited[i][x][y][z]) {
                                return 0;
                            }
                        }
                    }
                }
            }
            if(f != inputs.f[i]) {
                return 0;
            }
            if(r != inputs.r[i]) {
                return 0;
            }
        }

        double sum = 0.0;
        for(int i = 0; i < n; i++) {
            if(pos[0][i].size() == 0 && pos[1][i].size() == 0) {
                return 0;
            } else if(pos[0][i].size() == 0 || pos[1][i].size() == 0) {
                sum += (pos[0][i].size() + pos[1][i].size());
            } else if(is_same(pos[0][i], pos[1][i])) {
                sum += 1.0 / pos[0][i].size();
            } else {
                return 0;
            }
        }
        long long score = round(1e9 * sum); 
        return score;
    }

};

namespace Solver {

    vector<vector<tuple<int,int,int>>> generateBaseModels() {
        vector<vector<tuple<int,int,int>>> res(2);
        for(int i = 0; i < 2; i++) {
            for(int x = 0; x < inputs.D; x++) {
                for(int y = 0; y < inputs.D; y++) {
                    for(int z = 0; z < inputs.D; z++) {
                        if(inputs.f[i][z][x] == 1 && inputs.r[i][z][y] == 1) {
                            res[i].push_back({x, y, z});
                        }
                    }
                }
            }
        }
        return res;
    }

    Output testSolve() {
        auto baseModels = Solver::generateBaseModels();
        Output outputs;
        outputs.n = 0;
        for(int i = 0; i < 2; i++) {
            for(auto [x, y, z]: baseModels[i]) {
                outputs.ans[i][x][y][z] = ++outputs.n;
            }
        }
        return outputs;
    }

    void generateInitialPos(int initial_points, int b, vector<tuple<int,int,int>>& initial_pos) {

        set<tuple<int,int,int>> seen;
        for(auto e : initial_pos) seen.insert(e);

        if(initial_pos.size() == 0) {
            auto v = out8[b];
            shuffle(v.begin(), v.end(), Ryuka.engine);
            for(int p = 1; p <= min(8, initial_points); p++) {
                if(seen.count(v[p-1])) continue;
                initial_pos.push_back(v[p-1]);
                seen.insert(v[p-1]);
            }
        }


        int pp = initial_pos.size() + 1;

        #if 0
        vector<tuple<int,int,int>> cands;
        int pos = 0;
        for(int p = pp; p <= initial_points; p++) {
            while(true) {
                auto [x, y, z] = outercells[b][pos];
                if(!seen.count({x, y, z})) {
                    cands.push_back({x, y, z});
                    seen.insert({x, y, z});
                    break;
                } else {
                    pos += 1;
                    pos %= outercells[b].size();
                }
            }
        }
        shuffle(cands.begin(), cands.end(), Ryuka.engine);
        for(int p = pp; p <= initial_points; p++) {
            initial_pos.push_back(cands[p-pp]);  
        }
        #endif

        for(int p = pp; p <= initial_points; p++) {
            while(true) {
                int x = Ryuka.rand(inputs.D);
                int y = Ryuka.rand(inputs.D);
                int z = Ryuka.rand(inputs.D);
                if(vaild_pos[b][x][y][z] && !seen.count({x, y, z})) {
                    initial_pos.push_back({x, y, z});
                    seen.insert({x, y, z});
                    break;
                }
            }
        }
    }

    Output randomBFSSolve(int initial_points, int base, vector<vector<tuple<int,int,int>>>& initial_pos) {

        /*変数宣言*/
        int rbase = 1 - base;
        vector assigned(2, vector(inputs.D, vector(inputs.D, vector<int>(inputs.D, -1))));
        vector<vector<tuple<int,int,int>>> shape_b(initial_points+1), shape_rb(initial_points+1);
        vector<vector<tuple<int,int,int>>> cands_b(initial_points+1), cands_rb(initial_points+1);        
        vector seen_cands(2, vector(initial_points+1, vector(inputs.D, vector(inputs.D, vector<int>(inputs.D, false)))));        
        auto add_cands = [&](vector<tuple<int,int,int>>& cands, int p, int b, int x, int y, int z) -> void {
            for(int k = 0; k < 6; k++) {
                auto [dx, dy, dz] = d3[k];
                int kx = x + dx;
                int ky = y + dy;
                int kz = z + dz;
                if(kx < 0 || ky < 0 || kz < 0 || kx >= inputs.D || ky >= inputs.D || kz >= inputs.D) continue;
                if(assigned[b][kx][ky][kz] == -1 && vaild_pos[b][kx][ky][kz] && !seen_cands[b][p][kx][ky][kz]) {
                    cands.push_back({kx, ky, kz});
                    seen_cands[b][p][kx][ky][kz] = true;
                }
            }
        };


        /*初期化*/
        queue<int> que;
        for(int p = 1; p <= initial_points; p++) {
            {
                que.push(p);
                auto [x, y, z] = initial_pos[base][p-1];
                assigned[base][x][y][z] = p;
                shape_b[p].push_back({x, y, z});
                add_cands(cands_b[p], p, base, x, y, z);
            }
            {
                auto [x, y, z] = initial_pos[rbase][p-1];
                assigned[rbase][x][y][z] = p;
                shape_rb[p].push_back({x, y, z});
                add_cands(cands_rb[p], p, rbase, x, y, z);
            }
        }

        /* BFS */
        while(!que.empty()) {
            int p = que.front();
            que.pop();
            shuffle(cands_b[p].begin(), cands_b[p].end(), Ryuka.engine);
            for(auto [bx, by, bz] : cands_b[p]) {
                if(assigned[base][bx][by][bz] != -1) continue;
                bool f = false;
                auto tmp_shape_b = shape_b[p];
                tmp_shape_b.push_back({bx, by, bz});
                shuffle(cands_rb[p].begin(), cands_rb[p].end(), Ryuka.engine);
                for(auto [rbx, rby, rbz] : cands_rb[p]) {
                    if(assigned[rbase][rbx][rby][rbz] != -1) continue;
                    auto tmp_shape_rb = shape_rb[p];
                    tmp_shape_rb.push_back({rbx, rby, rbz});
                    if(Utils::is_same(tmp_shape_b, tmp_shape_rb)) {
                        assigned[base][bx][by][bz] = p;
                        assigned[rbase][rbx][rby][rbz] = p;
                        que.push(p);
                        shape_b[p].push_back({bx, by, bz});
                        shape_rb[p].push_back({rbx, rby, rbz});
                        add_cands(cands_b[p], p, base, bx, by, bz);
                        add_cands(cands_rb[p], p, rbase, rbx, rby, rbz);
                        f = true;
                        break;
                    }
                }
                if(f) break;
            }
        }

        Output outputs;
        outputs.n = initial_points;
        outputs.ans = assigned;
        for(int i = 0; i < 2; i++) {
            for(int x = 0; x < inputs.D; x++) {
                for(int y = 0; y < inputs.D; y++) {
                    for(int z = 0; z < inputs.D; z++) {
                        if(outputs.ans[i][x][y][z] == -1) {
                            outputs.ans[i][x][y][z] = 0;
                        }
                    }
                }
            }
        }
        return outputs;
    }

    pair<Output, vector<vector<tuple<int,int,int>>>> fillBlanks(Output outputs, bool opt = true) {
        vector tmp_f(2, vector(inputs.D, vector<int>(inputs.D, false)));
        vector tmp_r(2, vector(inputs.D, vector<int>(inputs.D, false)));
        vector additional_pos(2, vector<tuple<int,int,int>>());
        for(int i = 0; i < 2; i++) {
            for(int x = 0; x < inputs.D; x++) {
                for(int y = 0; y < inputs.D; y++) {
                    for(int z = 0; z < inputs.D; z++) {
                        if(outputs.ans[i][x][y][z] > 0) {
                            tmp_f[i][z][x] = true;
                            tmp_r[i][z][y] = true;
                        }
                    }
                }
            }
        }
        vector<int> next_p(2, outputs.n+1);
        auto tmp_func = [&](bool ty) {
            for(int i = 0; i < 2; i++) {
                for(int x = 0; x < inputs.D; x++) {
                    for(int y = 0; y < inputs.D; y++) {
                        for(int z = 0; z < inputs.D; z++) {
                            bool f1 = (inputs.f[i][z][x] == 1 && inputs.r[i][z][y]);
                            bool f2 = outputs.ans[i][x][y][z] == 0;
                            bool f3 = !tmp_f[i][z][x];
                            bool f4 = !tmp_r[i][z][y];
                            bool f5 = (ty ? (f3 && f4) : (f3 || f4));
                            if(f1 && f2 && f5) {
                                outputs.ans[i][x][y][z] = next_p[i];
                                next_p[i]++;
                                tmp_f[i][z][x] = true;
                                tmp_r[i][z][y] = true;
                                additional_pos[i].push_back({x, y, z});
                            }
                        }
                    }
                }
            }
        };
        tmp_func(true);
        tmp_func(false);

        if(opt && min(next_p[0], next_p[1])-1 >= outputs.n + 1) {
            int base = 0;
            int rbase = 1 - base;
            int pp = outputs.n + 1;
            int mn = min(next_p[0], next_p[1]) - 1;
            queue<int> que;
            auto assigned = outputs.ans;
            auto initial_pos = additional_pos;
            for(int i = 0; i < 2; i++) {
                for(int x = 0; x < inputs.D; x++) {
                    for(int y = 0; y < inputs.D; y++) {
                        for(int z = 0; z < inputs.D; z++) {
                            if(assigned[i][x][y][z] == 0 || assigned[i][x][y][z] > mn) {
                                assigned[i][x][y][z] = -1;
                            }
                        }
                    }
                }
            }
            vector<vector<tuple<int,int,int>>> shape_b(mn+1), shape_rb(mn+1);
            vector<vector<tuple<int,int,int>>> cands_b(mn+1), cands_rb(mn+1);        
            vector seen_cands(2, vector(mn+1, vector(inputs.D, vector(inputs.D, vector<int>(inputs.D, false)))));        
            auto add_cands = [&](vector<tuple<int,int,int>>& cands, int p, int b, int x, int y, int z) -> void {
                for(int k = 0; k < 6; k++) {
                    auto [dx, dy, dz] = d3[k];
                    int kx = x + dx;
                    int ky = y + dy;
                    int kz = z + dz;
                    if(kx < 0 || ky < 0 || kz < 0 || kx >= inputs.D || ky >= inputs.D || kz >= inputs.D) continue;
                    if(assigned[b][kx][ky][kz] == -1 && vaild_pos[b][kx][ky][kz] && !seen_cands[b][p][kx][ky][kz]) {
                        cands.push_back({kx, ky, kz});
                        seen_cands[b][p][kx][ky][kz] = true;
                    }
                }
            };
            for(int p = pp; p <= mn; p++) {
                {
                    que.push(p);
                    auto [x, y, z] = initial_pos[base][p-pp];
                    assigned[base][x][y][z] = p;
                    shape_b[p].push_back({x, y, z});
                    add_cands(cands_b[p], p, base, x, y, z);
                }
                {
                    auto [x, y, z] = initial_pos[rbase][p-pp];
                    assigned[rbase][x][y][z] = p;
                    shape_rb[p].push_back({x, y, z});
                    add_cands(cands_rb[p], p, rbase, x, y, z);
                }
            }
            /* BFS */
            while(!que.empty()) {
                int p = que.front();
//                assert(p <= mn);
                que.pop();
                shuffle(cands_b[p].begin(), cands_b[p].end(), Ryuka.engine);
                for(auto [bx, by, bz] : cands_b[p]) {
                    if(assigned[base][bx][by][bz] != -1) continue;
                    bool f = false;
                    auto tmp_shape_b = shape_b[p];
                    tmp_shape_b.push_back({bx, by, bz});
                    shuffle(cands_rb[p].begin(), cands_rb[p].end(), Ryuka.engine);
                    for(auto [rbx, rby, rbz] : cands_rb[p]) {
                        if(assigned[rbase][rbx][rby][rbz] != -1) continue;
                        auto tmp_shape_rb = shape_rb[p];
                        tmp_shape_rb.push_back({rbx, rby, rbz});
                        if(Utils::is_same(tmp_shape_b, tmp_shape_rb)) {
                            assigned[base][bx][by][bz] = p;
                            assigned[rbase][rbx][rby][rbz] = p;
                            que.push(p);
                            shape_b[p].push_back({bx, by, bz});
                            shape_rb[p].push_back({rbx, rby, rbz});
                            add_cands(cands_b[p], p, base, bx, by, bz);
                            add_cands(cands_rb[p], p, rbase, rbx, rby, rbz);
                            f = true;
                            break;
                        }
                    }
                    if(f) break;
                }
            }
            outputs.n = mn;
            outputs.ans = assigned;
            for(int i = 0; i < 2; i++) {
                for(int x = 0; x < inputs.D; x++) {
                    for(int y = 0; y < inputs.D; y++) {
                        for(int z = 0; z < inputs.D; z++) {
                            if(outputs.ans[i][x][y][z] == -1) {
                                outputs.ans[i][x][y][z] = 0;
                                //assert(outputs.ans[i][x][y][z] <= mn);
                            }
                        }
                    }
                }
            }
            auto res = Solver::fillBlanks(outputs, false);
            outputs = res.first;
            #if 0
            for(int i = 0; i < 2; i++) {
                for(int x = 0; x < inputs.D; x++) {
                    for(int y = 0; y < inputs.D; y++) {
                        for(int z = 0; z < inputs.D; z++) {
                            if(outputs.ans[i][x][y][z] > 0) {
                                tmp_f[i][z][x] = true;
                                tmp_r[i][z][y] = true;
                            }
                        }
                    }
                }
            }
            outputs.n = mn;
            outputs.write();
            next_p = vector<int>(2, outputs.n+1);
            vector tmp_f(2, vector(inputs.D, vector<int>(inputs.D, false)));
            vector tmp_r(2, vector(inputs.D, vector<int>(inputs.D, false)));
            auto tmp_func2 = [&](bool ty) {
                for(int i = 0; i < 2; i++) {
                    for(int x = 0; x < inputs.D; x++) {
                        for(int y = 0; y < inputs.D; y++) {
                            for(int z = 0; z < inputs.D; z++) {
                                bool f1 = (inputs.f[i][z][x] == 1 && inputs.r[i][z][y] == 1);
                                bool f2 = outputs.ans[i][x][y][z] == 0;
                                bool f3 = !tmp_f[i][z][x];
                                bool f4 = !tmp_r[i][z][y];
                                bool f5 = (ty ? (f3 && f4) : (f3 || f4));
                                if(f1 && f2 && f5) {
                                    outputs.ans[i][x][y][z] = next_p[i];
                                    next_p[i]++;
                                    tmp_f[i][z][x] = true;
                                    tmp_r[i][z][y] = true;
                                }
                            }
                        }
                    }
                }
            };
            tmp_func2(true);
            tmp_func2(false);
            if(min(next_p[0], next_p[1]) != outputs.n+1) {
                outputs.write();
                Sera.info(to_string(mn));
                Sera.info(to_string(next_p[0]));
                Sera.info(to_string(next_p[1]));
                Sera.info(to_string(outputs.n+1));
            }
            assert(min(next_p[0], next_p[1]) == outputs.n+1);
            #endif 
        }


        if(next_p[0] != next_p[1]) {
            int b = (next_p[0] < next_p[1] ? 1 : 0);
            int mn = next_p[1-b];
            atcoder::dsu uf(next_p[b]);
            for(int x = 0; x < inputs.D; x++) {
                for(int y = 0; y < inputs.D; y++) {
                    for(int z = 0; z < inputs.D; z++) {
                        if(outputs.ans[b][x][y][z] >= mn) {
                            for(int k = 0; k < 6; k++) {
                                auto [dx, dy, dz] = d3[k];
                                int kx = dx + x;
                                int ky = dy + y;
                                int kz = dz + z;
                                if(kx >= 0 && ky >= 0 && kz >= 0 && kx < inputs.D && ky < inputs.D && kz < inputs.D) {
                                    if(outputs.ans[b][kx][ky][kz] >= mn) {
                                        uf.merge(outputs.ans[b][kx][ky][kz], outputs.ans[b][x][y][z]);
                                    }
                                }
                            }
                        }
                    }
                }
            }
            int id = mn;
            vector<int> tb(next_p[b], -1);
            for(int x = 0; x < inputs.D; x++) {
                for(int y = 0; y < inputs.D; y++) {
                    for(int z = 0; z < inputs.D; z++) {
                        if(outputs.ans[b][x][y][z] >= mn) {
                            int ld = uf.leader(outputs.ans[b][x][y][z]);
                            if(tb[ld] == -1) {
                                tb[ld] = id;
                                id++;
                            }
                            outputs.ans[b][x][y][z] = tb[ld];
                        }
                    }
                }
            }
            next_p[b] = id;
        }
        outputs.n = max(next_p[0], next_p[1]) - 1;
        return {outputs, additional_pos};
    }

}

int main() {
    Toki.init();
    inputs.read();
    FunaQ.init();
    auto baseModels = Solver::generateBaseModels();
    vaild_pos.resize(2, vector(inputs.D, vector(inputs.D, vector<int>(inputs.D, false))));
    for(int i = 0; i < 2; i++) {
        for(auto [x, y, z]: baseModels[i]) {
            vaild_pos[i][x][y][z] = true;
        }
    }
    #if 0
    for(int i = 0; i < 2; i++) {
        vector<pair<int,tuple<int,int,int>>> tmp;
        for(int x = 0; x < inputs.D; x++) {
            for(int y = 0; y < inputs.D; y++) {
                for(int z = 0; z < inputs.D; z++) {
                    int count = 0;
                    for(int k = 0; k < 6; k++) {
                        auto [dx, dy, dz] = d3[k];
                        int kx = x + dx;
                        int ky = y + dy;
                        int kz = z + dz;
                        if(kx >= 0 && ky >= 0 && kz >= 0 && kx < inputs.D && ky < inputs.D && kz < inputs.D) {
                            if(!vaild_pos[i][x][y][z]) count++;
                        }
                    }
                    if(vaild_pos[i][x][y][z]) tmp.push_back({count, {x, y, z}});
                }
            }
        }
        sort(tmp.begin(), tmp.end(), [](auto a, auto b) { return a.first > b.first; });
        for(auto e : tmp) outercells[i].push_back(e.second);
    }
    Sera.info("outercells[0].size : " + to_string(outercells[0].size()));
    Sera.info("outercells[1].size : " + to_string(outercells[1].size()));
    #endif
    vector<tuple<int,int,int>> B8 =  {{0, 0, 0}, {inputs.D, inputs.D, inputs.D}, {0, 0, inputs.D}, {inputs.D, inputs.D, 0},
         {0, inputs.D, 0}, {inputs.D, 0, inputs.D}, {inputs.D, 0, 0}, {0, inputs.D, inputs.D}};
    vector<vector<int>> min_dists(2, vector<int>(8, 1<<30));
    for(int i = 0; i < 2; i++) {
        for(int x = 0; x < inputs.D; x++) {
            for(int y = 0; y < inputs.D; y++) {
                for(int z = 0; z < inputs.D; z++) {
                    if(vaild_pos[i][x][y][z]) {
                        for(int k = 0; k < 8; k++) {
                            auto [kx, ky, kz] = B8[k];
                            int dx = abs(x - kx);
                            int dy = abs(y - ky);
                            int dz = abs(z - kz);
                            int dist = dx + dy + dz;
                            if(dist < min_dists[i][k]) {
                                min_dists[i][k] = dist;
                                out8[i][k] = {x, y, z};
                            }
                        }
                    }
                }
            }
        }
    }
    Record best;
    best.score = 1LL<<60;
    int num_iter = 0, num_swap = 0;
    float ave_cost = 0;
    vector<vector<tuple<int,int,int>>> prev_pos;
    int npool = 2*inputs.D;
    bool fpool = true;
    vector<Record> rec_pool;
    int ini_pts = 0;
    int uprange = 2*inputs.D;
    while(Toki.elapsed() + ave_cost < 5.3) {
        num_iter++;
        Record cur;
        if(num_iter < 3*npool) {
            cur.base = Ryuka.rand(2);
            cur.initial_points = inputs.D - 2 + (ini_pts++) % uprange;
            if(cur.initial_points < 3) cur.initial_points = 3;
            if(cur.initial_points >= 2*inputs.D) cur.initial_points = 2*inputs.D;
            Solver::generateInitialPos(cur.initial_points, cur.base, cur.initial_pos[cur.base]);
            Solver::generateInitialPos(cur.initial_points, 1 - cur.base, cur.initial_pos[1 - cur.base]);
        } else {
            if(rec_pool.size() >= npool+1) {
                sort(rec_pool.begin(), rec_pool.begin()+npool+1, [](auto a, auto b){return a.score < b.score;});
            }
            int cur_id = Ryuka.rand(min<int>(npool, rec_pool.size()));
            cur = rec_pool[cur_id];
            int proc = Ryuka.rand(3);
            for(int i = 0; i < 2; i++) {
                if(proc == 0) {
                    shuffle(cur.initial_pos[i].begin(), cur.initial_pos[i].end(), Ryuka.engine);
                    cur.initial_pos[i].pop_back();
                    Solver::generateInitialPos(cur.initial_points, i, cur.initial_pos[i]);
                } else if(proc == 1) {
                    if(i == 0) cur.initial_points++;
                    Solver::generateInitialPos(cur.initial_points, i, cur.initial_pos[i]);
                } else {
                    shuffle(cur.initial_pos[i].begin(), cur.initial_pos[i].end(), Ryuka.engine);
                    cur.initial_pos[i].pop_back();
                    if(i == 0) cur.initial_points--;
                    Solver::generateInitialPos(cur.initial_points, i, cur.initial_pos[i]);
                }
            }
        }
        auto outputs_ = Solver::randomBFSSolve(cur.initial_points, cur.base, cur.initial_pos);
        auto res = Solver::fillBlanks(outputs_);
        cur.outputs = res.first;
        cur.score = Utils::compute_score(cur.outputs);
        rec_pool.push_back(cur);
        if(rec_pool.size() >= npool+1) swap(rec_pool[rec_pool.size()-1], rec_pool[npool]);
        if(rec_pool.size() >= 100) {
            sort(rec_pool.begin(), rec_pool.end(), [](Record a, Record b){ return a.score < b.score; });
            auto tmp = rec_pool;
            rec_pool.resize(npool);
            for(int i = 0; i < npool; i++) rec_pool[i] = tmp[i];
        }
        if(cur.score < best.score) {
            num_swap++;
            best = cur;
        }
        ave_cost = Toki.elapsed() / num_iter;
    }
    best.outputs.write();
}   
