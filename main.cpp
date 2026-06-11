#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>
#include <numeric>
#include <cmath>
#include <iomanip>
#include <climits>
#include <string>
using namespace std;

mt19937 rng(chrono::steady_clock::now().time_since_epoch().count());

// Solucion: permutacion + makespan
struct Solution { vector<int> perm; int makespan; };

// Tiempos detallados de la programacion PFSP
struct ResultPFSP { int makespan; vector<vector<int>> inicio, fin; };

// C[i][j] = max(C[i-1][j], C[i][j-1]) + p[perm[i]][j]
ResultPFSP calcularMakespan(const vector<int>& perm, const vector<vector<int>>& p, int n, int m) {
    vector<vector<int>> C(n, vector<int>(m,0)), inicio(n, vector<int>(m,0)), fin(n, vector<int>(m,0));
    for (int i = 0; i < n; i++)
        for (int j = 0; j < m; j++) {
            inicio[i][j] = max(i>0 ? C[i-1][j] : 0, j>0 ? C[i][j-1] : 0);
            fin[i][j] = C[i][j] = inicio[i][j] + p[perm[i]][j];
        }
    return {C[n-1][m-1], inicio, fin};
}

// Heuristica NEH: ordena por tiempo total desc. e inserta minimizando makespan
Solution generarSolucionNEH(const vector<vector<int>>& p, int n, int m) {
    vector<int> orden(n);
    iota(orden.begin(), orden.end(), 0);
    sort(orden.begin(), orden.end(), [&](int a, int b) {
        int sa=0, sb=0;
        for (int j=0; j<m; j++) { sa+=p[a][j]; sb+=p[b][j]; }
        return sa > sb;
    });
    vector<int> parcial;
    for (int k = 0; k < n; k++) {
        int trabajo = orden[k], mejorMakespan = INT_MAX, mejorPos = 0;
        for (int pos = 0; pos <= (int)parcial.size(); pos++) {
            vector<int> tmp = parcial;
            tmp.insert(tmp.begin()+pos, trabajo);
            int ms = calcularMakespan(tmp, p, (int)tmp.size(), m).makespan;
            if (ms < mejorMakespan) { mejorMakespan = ms; mejorPos = pos; }
        }
        parcial.insert(parcial.begin()+mejorPos, trabajo);
    }
    Solution s; s.perm = parcial;
    s.makespan = calcularMakespan(s.perm, p, n, m).makespan;
    return s;
}

// Permutacion aleatoria para diversidad inicial
Solution permutacionAleatoria(const vector<vector<int>>& p, int n, int m) {
    Solution s; s.perm.resize(n);
    iota(s.perm.begin(), s.perm.end(), 0);
    shuffle(s.perm.begin(), s.perm.end(), rng);
    s.makespan = calcularMakespan(s.perm, p, n, m).makespan;
    return s;
}

// Swap: intercambia dos posiciones aleatorias
Solution aplicarSwap(const Solution& s, int n) {
    Solution ns = s;
    uniform_int_distribution<int> dist(0, n-1);
    int i = dist(rng), j = dist(rng);
    while (j==i) j = dist(rng);
    swap(ns.perm[i], ns.perm[j]);
    return ns;
}

// Insert: extrae trabajo en i y lo reinserta en j
Solution aplicarInsert(const Solution& s, int n) {
    Solution ns = s;
    uniform_int_distribution<int> dist(0, n-1);
    int i = dist(rng), j = dist(rng);
    while (j==i) j = dist(rng);
    int trabajo = ns.perm[i];
    ns.perm.erase(ns.perm.begin()+i);
    if (j>i) j--;
    ns.perm.insert(ns.perm.begin()+j, trabajo);
    return ns;
}

// Reverse: invierte subsegmento [i,j] (2-opt)
Solution aplicarReverse(const Solution& s, int n) {
    Solution ns = s;
    uniform_int_distribution<int> dist(0, n-1);
    int i = dist(rng), j = dist(rng);
    if (i>j) swap(i,j);
    reverse(ns.perm.begin()+i, ns.perm.begin()+j+1);
    return ns;
}

Solution aplicarOperador(int op, const Solution& s, int n) {
    if (op==0) return aplicarSwap(s, n);
    if (op==1) return aplicarInsert(s, n);
    return aplicarReverse(s, n);
}

// Parametros del Fireworks Algorithm
struct Parametros { int numFireworks, maxSparks, minSparks, gaussSparks, MAX_ITER; double ampMax; };

// FWA para PFSP: fw[0]=NEH, resto aleatorios; itera generando chispas
Solution fireworksAlgorithm(int n, int m, const vector<vector<int>>& p,
                             const Parametros& params, vector<int>& convergencia,
                             const string& nombreInstancia) {
    vector<Solution> fuegos;
    fuegos.push_back(generarSolucionNEH(p, n, m));
    for (int i=1; i<params.numFireworks; i++)
        fuegos.push_back(permutacionAleatoria(p, n, m));

    Solution mejor = *min_element(fuegos.begin(), fuegos.end(),
        [](const Solution& a, const Solution& b){ return a.makespan < b.makespan; });
    convergencia.clear();
    convergencia.push_back(mejor.makespan);

    uniform_int_distribution<int> opDist(0, 2);
    for (int iter=0; iter<params.MAX_ITER; iter++) {
        int msMin = fuegos[0].makespan, msMax = fuegos[0].makespan;
        for (auto& fw : fuegos) { msMin=min(msMin,fw.makespan); msMax=max(msMax,fw.makespan); }

        vector<Solution> chispas;
        for (auto& fw : fuegos) {
            double calidad = (msMax==msMin) ? 1.0 :
                (double)(msMax-fw.makespan+1e-6)/(double)(msMax-msMin+1e-6);
            int numChispas = max(params.minSparks, min(params.maxSparks,
                (int)(params.minSparks+(params.maxSparks-params.minSparks)*calidad)));
            int amplitud = max(1, (int)(params.ampMax*n*(1.0-calidad+0.1)));

            for (int s=0; s<numChispas; s++) {
                Solution chispa = fw;
                for (int a=0; a<amplitud; a++) chispa = aplicarOperador(opDist(rng), chispa, n);
                chispa.makespan = calcularMakespan(chispa.perm, p, n, m).makespan;
                chispas.push_back(chispa);
            }
            normal_distribution<double> gauss(0.0, 1.0);
            for (int g=0; g<params.gaussSparks; g++) {
                Solution chispa = fw;
                int nOps = max(1, (int)abs(gauss(rng)*n/4.0));
                for (int a=0; a<nOps; a++) chispa = aplicarSwap(chispa, n);
                chispa.makespan = calcularMakespan(chispa.perm, p, n, m).makespan;
                chispas.push_back(chispa);
            }
        }

        vector<Solution> pool;
        for (auto& fw : fuegos)  pool.push_back(fw);
        for (auto& ch : chispas) pool.push_back(ch);
        sort(pool.begin(), pool.end(),
             [](const Solution& a, const Solution& b){ return a.makespan < b.makespan; });

        vector<Solution> siguiente; siguiente.push_back(pool[0]);
        uniform_int_distribution<int> idxDist(1, (int)pool.size()-1);
        while ((int)siguiente.size() < params.numFireworks)
            siguiente.push_back(pool[idxDist(rng)]);
        fuegos = siguiente;

        Solution iterBest = *min_element(fuegos.begin(), fuegos.end(),
            [](const Solution& a, const Solution& b){ return a.makespan < b.makespan; });
        if (iterBest.makespan < mejor.makespan) mejor = iterBest;
        convergencia.push_back(mejor.makespan);
    }
    return mejor;
}

// Carga instancia desde archivo: formato A (m valores/fila) o B (pares col-val)
bool cargarInstancia(const string& archivo, int& n, int& m, vector<vector<int>>& p) {
    ifstream fin(archivo);
    if (!fin) { cerr << "No se pudo abrir el archivo: " << archivo << "\n"; return false; }
    fin >> n >> m;
    p.assign(n, vector<int>(m, 0));
    for (int i=0; i<n; i++) {
        string linea;
        while (getline(fin, linea))
            if (!linea.empty() && linea.find_first_not_of(" \t\r\n") != string::npos) break;
        istringstream ss(linea);
        vector<int> tokens; int val;
        while (ss >> val) tokens.push_back(val);
        if ((int)tokens.size() == m) {
            for (int j=0; j<m; j++) p[i][j] = tokens[j];
        } else if ((int)tokens.size() == 2*m) {
            for (int k=0; k<2*m; k+=2) {
                int col=tokens[k], tiempo=tokens[k+1];
                if (col>=0 && col<m) p[i][col] = tiempo;
            }
        } else {
            cerr << "Formato inesperado en fila " << i << ": " << tokens.size()
                 << " tokens (esperados " << m << " o " << 2*m << ")\n";
            return false;
        }
    }
    fin.close(); return true;
}

// Procesa una instancia: 10 corridas + tabla resumen + exporta convergencia y Gantt
void procesarInstancia(const string& archivo, const string& nombreInstancia) {
    int n, m; vector<vector<int>> p;
    cout << "\n" << string(60,'=') << "\n";
    cout << "INSTANCIA: " << nombreInstancia << "  (" << archivo << ")\n";
    cout << string(60,'=') << "\n";
    if (!cargarInstancia(archivo, n, m, p)) {
        cout << "[OMITIDA] No se encontro el archivo: " << archivo << "\n"; return;
    }
    cout << "Trabajos: " << n << "  |  Maquinas: " << m << "\n";

    Parametros params;

    params.numFireworks = 10;
    params.maxSparks    = 50;
    params.minSparks    = 2;
    params.ampMax       = 1.0;
    params.gaussSparks  = 5;
    params.MAX_ITER     = 0;
    if      (n<=6)  params.MAX_ITER = 300;
    else if (n<=12) params.MAX_ITER = 500;
    else            params.MAX_ITER = 1000;

    cout << "\nParametros utilizados:\n"
         << "  Fuegos artificiales : " << params.numFireworks << "\n"
         << "  Chispas min/max     : " << params.minSparks << " / " << params.maxSparks << "\n"
         << "  Chispas gaussianas  : " << params.gaussSparks << "\n"
         << "  Amplitud maxima     : " << params.ampMax << "\n"
         << "  Iteraciones maximas : " << params.MAX_ITER << "\n";

    int runs = 10;
    vector<int> makespans; vector<double> tiempos;
    Solution mejorGlobal; mejorGlobal.makespan = INT_MAX;
    vector<int> mejorConvergencia;

    cout << "\n" << string(60,'-') << "\n"
         << left << setw(8)<<"Corrida" << setw(10)<<"Makespan"
         << setw(12)<<"Tiempo (s)" << setw(12)<<"Iteraciones" << "Mejor secuencia\n"
         << string(60,'-') << "\n";

    for (int r=0; r<runs; r++) {
        auto t0 = chrono::high_resolution_clock::now();
        vector<int> convergencia;
        Solution sol = fireworksAlgorithm(n, m, p, params, convergencia, nombreInstancia);
        if (sol.makespan < mejorGlobal.makespan) { mejorGlobal=sol; mejorConvergencia=convergencia; }
        double seg = chrono::duration<double>(chrono::high_resolution_clock::now()-t0).count();
        makespans.push_back(sol.makespan); tiempos.push_back(seg);

        string seq = "[";
        for (int i=0; i<n; i++) { seq+="J"+to_string(sol.perm[i]+1); if(i<n-1) seq+=","; }
        seq += "]";
        cout << left << setw(8)<<(r+1) << setw(10)<<sol.makespan
             << setw(12)<<fixed<<setprecision(3)<<seg << setw(12)<<params.MAX_ITER << seq << "\n";
    }

    int mejorMs  = *min_element(makespans.begin(), makespans.end());
    int peorMs   = *max_element(makespans.begin(), makespans.end());
    double prom  = accumulate(makespans.begin(), makespans.end(), 0.0) / runs;
    double var   = 0; for (int v : makespans) var += (v-prom)*(v-prom);
    double desv  = sqrt(var/(runs-1));
    double tProm = accumulate(tiempos.begin(), tiempos.end(), 0.0) / runs;

    string mejorSeq = "[";
    for (int i=0; i<n; i++) { mejorSeq+="J"+to_string(mejorGlobal.perm[i]+1); if(i<n-1) mejorSeq+=","; }
    mejorSeq += "]";

    cout << "\n--- Tabla resumen (" << runs << " corridas) ---\n"
         << "Mejor makespan    : " << mejorMs << "\n"
         << "Peor makespan     : " << peorMs  << "\n"
         << "Promedio          : " << fixed<<setprecision(2)<<prom << "\n"
         << "Desv. estandar    : " << fixed<<setprecision(2)<<desv << "\n"
         << "Tiempo promedio   : " << fixed<<setprecision(3)<<tProm << "s\n"
         << "Mejor secuencia   : " << mejorSeq << "\n";

    string nombreArchivoConv = "convergencia_" + nombreInstancia + ".csv";
    ofstream conv(nombreArchivoConv);
    if (conv.is_open()) {
        conv << "Iteracion,Makespan\n";
        for (size_t i=0; i<mejorConvergencia.size(); i++)
            conv << i+1 << "," << mejorConvergencia[i] << "\n";
        conv.close();
        cout << ">> Archivo de convergencia exportado: " << nombreArchivoConv << "\n";
    }

    ResultPFSP datosGantt = calcularMakespan(mejorGlobal.perm, p, n, m);
    string nombreArchivoGantt = "gantt_" + nombreInstancia + ".csv";
    ofstream gantt(nombreArchivoGantt);
    if (gantt.is_open()) {
        gantt << "trabajo,maquina,inicio,fin,duracion\n";
        for (int i=0; i<n; i++)
            for (int j=0; j<m; j++)
                gantt << "J"<<(mejorGlobal.perm[i]+1) << ",M"<<(j+1) << ","
                      << datosGantt.inicio[i][j] << "," << datosGantt.fin[i][j] << ","
                      << (datosGantt.fin[i][j]-datosGantt.inicio[i][j]) << "\n";
        gantt.close();
        cout << ">> Archivo Gantt exportado con exito: " << nombreArchivoGantt << "\n";
    }
    cout << "Makespan final del Gantt: " << datosGantt.makespan << "\n";
}

// Main: procesa automaticamente las 3 instancias del proyecto
int main() {
    cout << "=== Fireworks Algorithm - PFSP ===\n"
         << "Grupo 4 | Algoritmo de Fuegos Artificiales\n"
         << "Procesando las 3 instancias...\n";

    vector<pair<string,string>> instancias = {
        {"instancia1_bas1.txt",  "Pequena_5x4" },
        {"instancia2_car5.txt",  "Mediana_10x6"},
        {"instancia3_reC01.txt", "Grande_20x5" }
    };
    for (auto& inst : instancias)
        procesarInstancia(inst.first, inst.second);

    cout << "\n" << string(60,'=') << "\nEjecucion completada.\n";
    return 0;
}
