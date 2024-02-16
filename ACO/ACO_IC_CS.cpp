#include<bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef unsigned long long ull;
#define forinc(i,a,b) for(int i=a,_b=b;i<=_b;++i)
#define fordec(i,a,b) for(int i=a;i>=b;--i)
#define forv(a,b) for(auto &a:b)
#define timer 1.0*clock()/CLOCKS_PER_SEC
#define pb push_back
#define double long double
#define all(a) a.begin(),a.end()
const int loop = 100, ant = 10, T = 20;
const int alp_size = 4,MAX = 1e9;
int num[1000],co[4],n,w;
int length;
double best_inf,ESP = 1e-9;
double rho = 0.5, t_max = 1, t_min = 1.0/1000;
vector<string> seq;
vector<int> best_pos;
vector<double> background;
vector<vector<int> > occurrences,b_o;
vector<vector<double> > pheromone_vertex;
vector<vector<double> > heuristic;
vector<vector<int> > sum_co;
vector<int> num_pos;
vector<vector<vector<int> > > f;
ll Rand(ll l,ll r) {
    return l + ((ll)rand() * (RAND_MAX + 1) * (RAND_MAX + 1) * (RAND_MAX + 1) +
                (ll)rand() * (RAND_MAX + 1) * (RAND_MAX + 1) +
                (ll)rand() * (RAND_MAX + 1) +
                rand())%(r-l+1);
}
double complex_score(vector<vector<int>> occurrences) {
        double score = 0;
        vector<int> e(alp_size);
        forinc(i, 0, w - 1) {
            int j = 0;
            forinc(k, 1, 3) if(occurrences[k][i] > occurrences[j][i]) j = k;
            e[j]++;
        }
        forinc(i, 0, 3) if(e[i]) {
            double p = (double)e[i]/w;
            score -= p*log2(p);
        }
        return score;
}
double information_content(int &n, // number of motif instances
                           const int &alp_size,
                           const int &w,
                           const vector<vector<int>> &occurrences,
                           const vector<double> &background) {
    double pseudo_count = 0;
    double result = 0.0;
    for (int i = 0; i < alp_size; ++i)
        for (int j = 0; j < w; ++j) if(occurrences[i][j]) {
            double frequency = (double) (occurrences[i][j] + pseudo_count) / (n + pseudo_count*alp_size);
            result += frequency*log2(frequency/background[i]);
        }
    return result;
}
double get_inf() {
    double IC = information_content( n, alp_size, w, occurrences, background );
    double CS = complex_score(occurrences);
    //cerr << IC <<" "<<CS <<"\n";
    double p = 0.8;
    return IC/w * p + CS * (1.0 - p);
}
void Add(int i,int j,vector<int> &pos) {
    ++n;
    pos[i] = j;
    forinc(p, j, j + w -1) occurrences[num[seq[i][p]]][p-j]++;
}
void Rev(int i,int j,vector<int> &pos) {
    --n;
    forinc(p, j, j + w -1) occurrences[num[seq[i][p]]][p-j]--;
}
void Change(int i, int j, int k,vector<int> &pos) {
    Rev(i, j, pos);
    Add(i, k, pos);
}
double get_IC(vector<vector<int> > &pos) {
    n = 0;
    vector<int> p(seq.size());
    forinc(i, 0, alp_size - 1) {
        forinc(j, 0, w - 1) occurrences[i][j] = 0;
    }
    forinc(i, 0, seq.size() - 1) {
        forinc(j, 0, length - 1) if(pos[i][j]) Add(i,j,p);
    }
    return get_inf();
}
bool cmp(vector<vector<int> > a,vector<vector<int> > b) {
    return get_IC(a) > get_IC(b);
}
void Enter() {
    memset( num, -1, sizeof(num));
    num['a'] = num['A'] = num['n'] = 0;
    num['c'] = num['C'] = 1;
    num['g'] = num['G'] = 2;
    num['t'] = num['T'] = 3;
    string s,S;
    cin >> w;
    while(cin >> s) {
        if(s[0] == '>')
        {
            if( !S.empty()) {
                seq.pb( S);
                S.clear();
            }
            continue;
        }
        else if( num[s[0]] >= 0) S += s;
    }
    seq.pb(S);
    length = seq[0].size();
    forinc(i, 0, alp_size - 1) {
        vector<int> occur( w );
        occurrences.pb( occur );
    }
    forinc(i, 0, seq.size() - 1) {
        vector<double> vertex( length - w +1, 1);
        vector<double> heuris( length - w +1);
        pheromone_vertex.pb( vertex );
        heuristic.pb(heuris);
    }
    forinc(i, 0, alp_size - 1) background.pb( 0.25 );
    return;
    int total = 0;
    forv(i, seq) {
        forv(j, i) co[num[j]]++;
        total += i.size();
    }
    forinc(i, 0, alp_size - 1)
    {
        background[i] = (double)co[i]/total;
        //cerr << background[i] << "\n";
    }
}
pair< int, int > next_random(pair< int, int > here,vector<int> &pos) {
    double total = 0, t_here = 0;
    forinc(i, 0, length - w) if(seq[here.first + 1][i]!='n') {
        Add(here.first + 1 , i, pos);
        heuristic[here.first + 1][i] = get_inf();
        total += heuristic[here.first + 1][i] *  pheromone_vertex[here.first + 1][i];
        Rev(here.first + 1, i, pos);
    }
    double random = (double) Rand(0,MAX)/MAX;
    forinc(i, 0, length - w) if(seq[here.first + 1][i]!='n') {
        t_here += heuristic[here.first + 1][i] *  pheromone_vertex[here.first + 1][i];
        if( random < t_here/total ) return {here.first + 1, i};
    }
}
vector<int> go_random() {
    vector<int> pos;
    pos.resize(seq.size());
    forinc(i, 0, alp_size - 1) {
        forinc(j, 0, w - 1) occurrences[i][j] = 0;
    }
    n = 0;
    pair< int, int > here = {-1, 0};
    while(here.first < (int)seq.size() - 1) {
        here = next_random( here, pos );
        Add(here.first, here.second, pos);
    }
    return pos;
}
void update_pheromone_vertex (vector<int> &pos) {
	forinc(i, 0, pos.size() - 1) {
	    int j = pos[i];
	    pheromone_vertex[i][j] = (1 - rho) * pheromone_vertex[i][j] + rho * t_max;
	    double c = (double)(w-1)/w;
	    double t = t_max*c;
	    for(int k = j + 1; k <= length - w; ++k,t*=c) pheromone_vertex[i][k] = (1 - rho) * pheromone_vertex[i][k] + rho * t;
        t = t_max*(w-1)/(w);
	    for(int k = j - 1; k >= 0; --k,t*=c) pheromone_vertex[i][k] = (1 - rho) * pheromone_vertex[i][k] + rho * t;
		/*forinc(j, 0, length - w) {
		    if(pos[i] == j)  pheromone_vertex[i][j] = (1 - rho) * pheromone_vertex[i][j] + rho * t_max;
		    else  pheromone_vertex[i][j] = (1 - rho) * pheromone_vertex[i][j] + rho * t_min;
		}*/
	}
}
void local_search(vector<int> &pos) {
    bool stop = 0;
    while(!stop) {
        stop = 1;
        forinc(i, 0, seq.size() - 1) {
            forinc(j, 0, length - w) if(seq[i][j]!='n'){
                double p_inf = get_inf();
                int p = pos[i];
                Rev(i,p,pos),Add(i,j,pos);
                if(get_inf() > p_inf) stop = 0;
                else Rev(i,j,pos),Add(i,p,pos);
            }
        }
    }
}
vector<vector<int> > local_full(vector<int> &pos) {
    n = 0;
    forinc(i, 0, alp_size - 1) {
        forinc(j, 0, w - 1) occurrences[i][j] = 0;
    }
    vector<vector<int> > p;
    forinc(i, 0, pos.size() - 1) {
        vector<int> New(length);
        p.pb(New);
        p[i][pos[i]] = 1;
        Add(i,pos[i],pos);
    }
    return p;
    vector<int> block;
    while(1) {
        bool ok = 0;
        // Change
        forinc(i, 0, seq.size() - 1) {
            forinc(j, 0, length - w) if(!p[i][j]&&seq[i][j]!='n') {
                block.clear();
                forinc(k, max(0, j - w + 1), min(length - w, j + w -1)) if(p[i][k]) block.pb(k);
                if(block.size()) {
                    double pre_inf = get_inf();
                    forv(k,block) p[i][k] = 0,Rev(i,k, pos);
                    p[i][j] = 1,Add(i,j, pos);
                    if(get_inf() > pre_inf) ok=1;
                    else
                    {
                        forv(k,block) p[i][k] = 1,Add(i,k, pos);
                        p[i][j]=0,Rev(i,j, pos);
                    }
                }
                else if(block.size()==0) {
                    double pre_inf = get_inf();
                    forinc(k, 0, length - w) if(p[i][k]) {
                        p[i][k] = 0, p[i][j] = 1;
                        Change(i, k, j, pos);
                        if(get_inf() > pre_inf) {ok=1;break;}
                        else p[i][k] = 1, p[i][j] = 0,Change(i, j, k, pos);
                    }
                }
            }
        }
        // Add
        forinc(i, 0, seq.size() - 1) {
            forinc(j, 0, length - w) if(!p[i][j]&&seq[i][j]!='n') {
                bool block=0;
                forinc(k, max(0, j - w + 1), min(length - w, j + w -1)) if(p[i][k]) block=1;
                if(!block) {
                    double pre_inf = get_inf();
                    p[i][j] = 1;
                    Add(i, j, pos);
                    if(get_inf() > pre_inf) ok=1;
                    else p[i][j] = 0,Rev(i,j, pos);
                }
            }
        }
        // Rev
        forinc(i, 0, seq.size() - 1) {
            int res=0;
            forinc(j, 0, length - w) if(p[i][j]) res++;
            forinc(j, 0, length - w)  if(res>1&&p[i][j]) {
                double pre_inf = get_inf();
                p[i][j] = 0,Rev(i, j, pos);
                if(get_inf() > pre_inf) ok=1,--res;
                else p[i][j] = 1,Add(i,j, pos);
            }
        }
        if(!ok) break;
    }
    return p;
}
void Print() {
    cout << "========================================\n";
    forinc(i, 0, seq.size() - 1) {
        cout << i+1 << "=" <<best_pos[i] <<  "\n";
    }
    //cout << "time = " <<setprecision(6) << fixed << timer << "\n";
}
void Solve() {
    int best_not_increase = 0;
    forinc(i, 1, loop) {
        double l_inf = 0;
        vector<int> l_pos;
        vector<vector<int> > l_o;
        forinc(j ,1, ant) {
            vector<int> r_pos = go_random();
            local_search(r_pos);
            double r_inf = get_inf();
            if(r_inf > l_inf ) {
                l_inf = r_inf;
                l_pos = r_pos;
                //l_o = occurrences;
            }
            //f.pb(local_full(r_pos));
        }
        //cout << l_inf <<"\n";
        best_not_increase++;
        if(l_inf > best_inf) {
            best_inf = l_inf;
            best_pos = l_pos;
            best_not_increase = 0;
            //b_o = l_o;
        }
        if(best_not_increase == T) break;
        update_pheromone_vertex( best_pos );
        cerr << "loop " << i << " ended.\n";
        //Print();
        //exit(0);
    }
}
void Print_best() {
    vector<vector<int> > pos = local_full(best_pos);
    cout << "IC = "<<get_IC(pos) << "\n";
    forinc(i, 0, pos.size() - 1) {
        cout << i+1 << "=";
        forinc(j, 0, length - w) if(pos[i][j]) cout << j << " ";
        cout << "\n";
    }
    cout << "time = " <<setprecision(6) << fixed << timer << "\n";
}
int main() {
    ios_base::sync_with_stdio(0);
    cin.tie(0);
    freopen("test.inp","r",stdin);
    freopen("test.out","w",stdout);
    srand(time(NULL));
    Enter();
    Solve();
    Print_best();
    /*sort(all(f));
    f.erase(unique(all(f)),f.end());
    sort(all(f),cmp);
    int cnt = 0;
    cerr << "/////////////////////////////////////////////////\n" ;
    forv(pos, f) {
        cerr << ++cnt << "\n";
        cout << get_IC(pos) << "\n";
        forinc(i, 0, pos.size() - 1) {
            cout << i+1 << "=";
            forinc(j, 0, length - w) if(pos[i][j]) cout << j << " ";
             cout << "\n";
        }
        cerr << "/////////////////////////////////////////////////\n" ;
    }*/
}

