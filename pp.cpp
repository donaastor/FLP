#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <vector>
#ifdef __linux__
#include "include/armadillo"
#else
#include <armadillo>
#endif
#undef max
#undef min

#ifndef __linux__
#pragma warning(disable:4996)
#endif
//using namespace std;
using namespace arma;
//nasumicno bira ceo broj iz intervala [0,i-1]
int randint(int i) {
	return randi<int>(distr_param(0, i - 1));
}
//nasumicno bira k elemenata iz skupa v
uvec randsample(uvec v, int k) {
	uvec d(k);
	if (k == 1) {
		int ind = randint(v.n_elem);
		d(0) = v(ind);
	}
	else {
		uvec v_shuffle = arma::shuffle(v);
		d = v_shuffle(span(0, k - 1));
	}
	return d;
}
//vraca vektor koji predstavlja razliku vektora v i u, tj. v\u
uvec setdiff(uvec v, urowvec u) {
	int l = v.n_elem;
	uvec s(l);
	int k = 0;
	for (int i = 0; i < l; i++) {
		if (all(u != v(i))) {
			s(k) = v(i);
			k += 1;
		}
	}
	if (k == 0) {
		return uvec();
	}
	return s(span(0, k - 1));
}
uvec setdiff2(uvec v, uvec u) {
	int l = v.n_elem;
	uvec s(l);
	int k = 0;
	for (int i = 0; i < l; i++) {
		if (all(u != v(i))) {
			s(k) = v(i);
			k += 1;
		}
	}
	if (k == 0) {
		return uvec();
	}
	return s(span(0, k - 1));
}
uvec tumbani(int n, int m, mat C) {
	uvec Cl = zeros<uvec>(m);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < m; j++) {
			if (C(i, j) != j + 1) {
				Cl(j) += 1;
			}
		}
	}
	return Cl;
}
/*
racuna funkciju cilja
res - resenje
C - matrica rastojanja
*/
int funkcija_cilja(uvec res, int n, int m, int p, mat C) {
	int f_res = 0;
	umat M = zeros<umat>(n, m);
	uvec indeksi = zeros<uvec>(p); //indeksi kolona od p, indeksi=[0 0]
	//uvec red = zeros<uvec>(p);
	//uvec minimumi = zeros<uvec>(n);
	//int mini = m; //min je najvise jednak broju kolona
	int l = 0;
	uvec red = zeros<uvec>(p); //red = [0 0]
	uvec suma = zeros<uvec>(p); //suma = [0 0]
	//int n=res.n_elem; //.n_elem - total number of elements
	for (int k = 0; k < m; k++) {
		if (res(k) == 1) { //res = [0 1 0 0 1]
			indeksi(l) = k; //indeksi=[1,4]
			l += 1;
		}
	}
	//int h = indeksi.n_elem;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < p; j++) { //C = [2 33 1 55 14]
			red(j) = C(i, indeksi(j)); //red = [33, 14]
		}
		uint32_t najmanji = red(0);
		uint32_t najmanji_indeks = 0;
		for (int k = 1; k < p; k++) {
			if (najmanji > red(k)) {
				najmanji = red(k);
				najmanji_indeks = k;
			}
		}
		M(i, indeksi(najmanji_indeks)) = 1;
	}
	//cout << "M: " << M << endl;
	for (int j = 0; j < p; j++)
		for (int i = 0; i < n; i++)
			suma(j) += M(i, indeksi(j));
	//cout << "maks " << max(suma) << endl;
	//cout << "min " << min(suma) << endl;
	f_res = max(suma) - min(suma);
	return f_res;
}
/*
azurira najbolje resenje ukoliko je doslo do poboljsanja
x_best - trenutno najbolje resenje
f_best - vr. funkcije cilja resenja x_best
F - vrednosti funkcija cilja svih (nekog broja) suseda koji se nalaze u okolini tekuceg resenja
P - susedi tekuceg resenja (po kolonama)
*/
void azuriraj_najbolje_resenje(uvec& y_best, int& f_best, vec F, umat P) {
	uword ind = index_min(F); //vraca indeks najmanje vrednosti
	uvec new_y_best = P.col(ind); //novi y koji je najbolji je u toj koloni
	int new_f_best = F(ind); //F od te kolone je nova najbolja funkcija cilja
	if (new_f_best < f_best) {
		f_best = new_f_best;
		y_best = new_y_best;
	}
}
/*
Za primenu operatora umetanja, zamene i inverzije potrebno je najpre
odrediti 2 pozicije u resenju, a potom primeni odgovarajuci operator
funkcija formira matricu I svih kombinacija indeksa (pozicija) pomocu kojih (uz primenu nekog od operatora)
mozemo generisati sve susede
*/
void sve_komb_ind(int m, umat& I) {
	int br = 0;
	uint32_t um = m;
	for (uint32_t i = 0; i < um; i++) {
		for (uint32_t j = 0; j < um; j++) {
			if (i != j) { //nista dijagonale
				uvec ind = { i,j };
				I.col(br) = ind;
				br += 1;
			}
		}
	}
}
/*
operator umetanja
y - resenje
ind - vektor duzine 2 (dva indeksa, ind(0) odredjuje trenutnu poziciju elementa koji treba postaviti (umetnuti) na poziciju ind(1))
*/
uvec insertion(uvec y, uvec ind) {
	int m = y.n_elem;
	//uvec v=linspace<uvec>(0,n-1,n);
	//uvec ind=randsample(v,2);
	//cout<<ind<<"****ind***"<<endl;
	//el. na poziciji ind(0) premestamo na poziciju ind(1)
	// int el = ind(0); //ind(0) = 4, ind(1) = 2      // ne koristi se
	uvec new_y = y;//new_y = y = [1 2 3 4 5]
	new_y.replace(y(ind(1)), y(ind(0))); //new_y = [1 2 5 4 3]
	if (ind(0) < ind(1)) {
		int t = ind(0) - 1;
		if (t < 0) {
			new_y = join_vert(y(span(ind(0) + 1, ind(1))), new_y(span(ind(1), m - 1)));
		}
		else {
			new_y = join_vert(join_vert(y(span(0, t)), y(span(ind(0) + 1, ind(1)))), new_y(span(ind(1), m - 1)));
		}
	}
	else {
		int s = ind(0) + 1; //s = 5
		if (s > m - 1) {
			new_y = join_vert(new_y(span(0, ind(1))), y(span(ind(1), ind(0) - 1)));
		} //new_y = [1 2 5 3 4]
		else {
			new_y = join_vert(join_vert(new_y(span(0, ind(1))), y(span(ind(1), ind(0) - 1))), y(span(s, m - 1)));
		}
	}
	return new_y;
}
/*
operator zamene
x - resenje
ind - vektor duzine 2 (elemnti na pozicijama ind(0) i ind(1) menjaju svoja mesta)
*/
uvec exchange(uvec y, uvec ind) {
	// int m = y.n_elem;      // ne koristi se
	//uvec v=linspace<uvec>(0,n-1,n);
	//uvec ind=randsample(v,2);
	//cout<<ind<<"****ind***"<<endl;
	//elementi na pozicijama ind(0) i ind(1) menjaju mesta
	int el0 = y(ind(0));
	int el1 = y(ind(1));
	uvec new_y = y;
	new_y(ind(0)) = el1;
	new_y(ind(1)) = el0;
	return new_y;
}
/*
operator inverzije
x - resenje
ind - vektor duzine 2 (deo resenja od pozicije ind(0) do ind(1) uzima se u obrnutom poretku)
*/
uvec inversion(uvec y, uvec ind) {
	int m = y.n_elem;
	//uvec v=linspace<uvec>(0,n-1,n);
	//uvec ind=randsample(v,2);
	//cout<<ind<<"****ind***"<<endl;
	int ind1 = min(ind);
	int ind2 = max(ind);
	int t = ind1 - 1;
	int s = ind2 + 1;
	uvec new_y;

	if (t<0 && s>(m - 1)) {
		new_y = reverse(y);
	}
	else if (t < 0) {
		new_y = join_vert(reverse(y(span(ind1, ind2))), y(span(s, m - 1)));
	}

	else if (s > (m - 1)) {
		new_y = join_vert(y(span(0, t)), reverse(y(span(ind1, ind2))));
	}
	else {
		new_y = join_vert(join_vert(y(span(0, t)), reverse(y(span(ind1, ind2)))), y(span(s, m - 1)));
	}
	return new_y;
}
/*Lokalna pretraga
y - tekuce resenje ciju okolinu pretrazujemo
f - vrednost funckije cilja resenja y
ind1 - indikator kojim se bira operator (0 - inverzija, 1 - zamena, 2 - umetanje)
ind2 -  indikator kojim se odredjuje startegija odabira suseda koji postaje novo tekuce resenje (0 - first improvement, 1 - best improvement)
C - matrica rastojanja
funkcija vraca lokalni optimum
*/
uvec local_search(uvec y, int n, int p, int f, int ind1, int ind2, mat C) {
	uvec y_best = y; //y_best = [0 1 0 0 1]
	int m = y.n_elem; //m = 5
	int f_best = f;
	umat I = zeros<umat>(2, m * m - m);
	sve_komb_ind(m, I); //kombinacija pozicija
	vec F = zeros<vec>(m * m - m);
	umat P = zeros<umat>(m, m * m - m);
	int indicator = 0;
	int indicator_FI = 0;
	uvec ind;
	int br;
	while (indicator == 0 && indicator_FI == 0) {
		if (ind2 == 0) {
			indicator_FI = 1;
		}
		for (br = 0; br < (m * m - m); br++) {
			ind = I.col(br);
			if (ind1 == 0) {
				P.col(br) = inversion(y_best, ind);
			}
			else if (ind1 == 1) {
				P.col(br) = exchange(y_best, ind);
			}
			else {
				P.col(br) = insertion(y_best, ind);
			}
			F(br) = funkcija_cilja(P.col(br), n, m, p, C);
			if (ind2 == 0 && F(br) < f_best) {
				//first improvement 
				indicator_FI = 0;
				azuriraj_najbolje_resenje(y_best, f_best, F(span(0, br)), P.cols(span(0, br)));
				break;
			}
		}
		if (ind2 == 1) {
			if (min(F) >= f_best) {
				indicator = 1;
				//cout<<"nije doslo do poboljsanja"<<endl;
			}
			else {
				azuriraj_najbolje_resenje(y_best, f_best, F, P);
			}
		}
	}
	return y_best;
}
/*
nasumicno-pohlepni algoritam za konstrukciju inicijalnog resenja GRASP metode
n - broj korisnika
m - broj trznih centara
p - broj otvenih trznih centara (u okviru m)
C - matrica rastojanja
funkcija vraca inicijalno resenje
*/
uvec konstrukcija_resenja_pohlepnim_alg(int n, int m, int p, mat C) {
	//uvec v=linspace<uvec>(0,n-1,n);
	uvec y = zeros<uvec>(m); //vektor y = [0 0 0 0 0], m = 5
	uvec v = linspace<uvec>(1, m, m);
	uvec new_el = randsample(v, 1);
	new_el(0) = new_el(0) - 1; //new_el = 1
	y(new_el(0)) = 1; //y = [0 1 0 0 0]
	uvec Cl;
	uvec Rcl;
	uvec indeksi;
	Cl = tumbani(n, m, C); //Cl = [0 3 0 1 2]
	Rcl = zeros<uvec>(m);
	for (int i = 1; i < p; i++) {
		//Cl=sort_index(C.col(x(i-1))); 
		//Cl=setdiff2(Cl,x(span(0,i-1)));
		Cl(new_el(0)) = -1; //Cl = [0 -1 0 1 2]
		int maks = max(Cl); //maks = 2
		Cl(new_el(0)) = m + 1; //Cl = [0 6 0 1 2]
		int mini = min(Cl); //min = 0
		int l = 0;

		for (int k = 0; k < m; k++) {
			if (Cl(k) >= mini + 0.5 * (maks - mini) && y(k) != 1) {
				Rcl(l) = k;
				l += 1;
			}
		}
		Rcl = Rcl(span(0, l - 1));
		new_el = randsample(Rcl, 1);
		y(new_el(0)) = 1;
	}
	return y;
}



using namespace std::chrono;
typedef steady_clock::time_point vreme;
inline vreme sada() {
	return steady_clock::now();
}
inline double merivreme(vreme poc, vreme kraj) {
	return duration_cast<duration<double>>(kraj - poc).count();
}

struct primer;
struct insta {
	primer* PR;
	uvec y;
	uint32_t f;
	uvec ly;
	uint32_t lf;
	int iter;
	int pit;
	uvec y_best;
	uint32_t f_best;
	uvec resenja;
	rowvec vremena;
	inline double ukvreme();
	inline double najvreme();
	inline void pisivreme(const char* txt);
	inline void pisivreme(const char* txt, int ix);
	inline void pocni(primer* group);
	inline void iteriraj();
};
struct primer {
	int n, m, p;
	int maxiter, it;
	mat C;
	constexpr static int INS = 15
		;	insta ie[INS];
	inline int ucitaj(const char* fname, double mik) {
		FILE* file = fopen(fname, "r");
		if (file == nullptr)
			return -1;
		fscanf(file, "%d", &m);
		fscanf(file, "%d", &n);
		fscanf(file, "%d", &p);
		maxiter = n * mik;
		C = zeros<mat>(n, m);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				int num;
				fscanf(file, "%d", &num);
				C(i, j) = num;
			}
		}
		fclose(file);
		return 0;
	}
	inline uint32_t maxres() {
		uint32_t najg = ie[0].f_best;
		for (int i = 1; i < INS; ++i)
			if (ie[i].f_best > najg)
				najg = ie[i].f_best;
		return najg;
	}
	inline uint32_t minres() {
		uint32_t najb = ie[0].f_best;
		for (int i = 1; i < INS; ++i)
			if (ie[i].f_best < najb)
				najb = ie[i].f_best;
		return najb;
	}
	inline void resi(int goal) {
		int ispiso = 0;
		static const char* tems = "                                                                                "; // 60
		static const char* zevs = "||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||"; // 60
		for (int i = 0; i < INS; ++i)
			ie[i].pocni(this);
		int res = maxres();
		if (res <= goal) {
			it = 0;
			printf("%s done\n", tems);
			return;
		}
		bool ng = false;
		int zavi = maxiter;
		for (int j = 1; j <= zavi; ++j) {
			for (int i = 0; i < INS; ++i)
				ie[i].iteriraj();
			int nisp = (80 * j) / maxiter;
			printf("%.*s", nisp - ispiso, zevs);
			fflush(stdout);
			ispiso = nisp;
			res = maxres();
			if ((ng == false) & (res <= goal)) {
				ng = true;
				zavi = j * 2;
				if (zavi > maxiter) zavi = maxiter;
				it = zavi;
			}
		}
		printf(" done\n");
	}
	inline double ukvreme() {
		double x = 0;
		for (int i = 0; i < INS; ++i)
			x += ie[i].ukvreme();
		return x;
	}
	inline double najvreme() {
		double x = 0;
		for (int i = 0; i < INS; ++i)
			x += ie[i].najvreme();
		return x;
	}
	inline int print(int ord, const char* name, int nlen, char* os) {
		int sh = 0;
		sh += sprintf(os + sh, "instance:   %2d, %.*s\r\nm,n,p:      %d %d %d\r\n", ord, nlen, name, m, n, p);
		sh += sprintf(os + sh, "\r\n       sol_i         t_i      ttot_i       gap_i    best_sol           t        ttot        agap       sigma        iter\r\n");
		int best_sol = minres();
		double avgt = najvreme(), avgttot = ukvreme(), ssol = 0, qsol = 0;
		for (int i = 0; i < INS; ++i)ssol += ie[i].f_best;
		for (int i = 0; i < INS; ++i)qsol += ie[i].f_best * ie[i].f_best;
		avgt /= INS; avgttot /= INS; ssol /= INS; qsol /= INS;
		const double kf = 100.0 / best_sol;
		qsol = kf * sqrt(qsol - ssol * ssol);
		ssol = kf * (ssol - best_sol);
		char gaps[13]; gaps[12] = 0;
		auto cgap = [&](int ind)->const char* {
			double g = kf * (ie[ind].f_best - best_sol);
			if (g != 0) sprintf(gaps, "%12.5f", g);
			else memcpy(gaps, "           -", 12);
			return gaps;
		};
		sh += sprintf(os + sh, "%12d%12.5f%12.5f%s%12d%12.5f%12.5f%12.5f%12.5f%12d\r\n",
			ie[0].f_best, ie[0].najvreme(), ie[0].ukvreme(), cgap(0),
			best_sol, avgt, avgttot, ssol, qsol, it);
		for (int i = 1; i < INS; ++i)
			sh += sprintf(os + sh, "%12d%12.5f%12.5f%s\r\n",
				ie[i].f_best, ie[i].najvreme(), ie[i].ukvreme(), cgap(i));
		sh += sprintf(os + sh, "------------------------------------------------------------------------------------------------------------------------\r\n");
		return sh;
	}
};
inline double insta::ukvreme() {
	return vremena(iter - 1);
}
inline double insta::najvreme() {
	return vremena(pit);
}
inline void insta::pisivreme(const char* txt, int ix) {
	char stocc[17] = "                ";
	int tl = strlen(txt);
	memcpy(stocc, txt, tl);
	stocc[tl] = ':';
	printf("%s%.2f seconds\n", stocc, vremena(ix));
}
inline void insta::pisivreme(const char* txt) {
	return pisivreme(txt, iter - 1);
}
inline void insta::pocni(primer* group) {
	PR = group;
	f = 0;
	resenja = zeros<uvec>(PR->maxiter + 1);
	vremena = zeros<rowvec>(PR->maxiter + 1);
	iter = 1;
	vreme t1 = sada();
	y = konstrukcija_resenja_pohlepnim_alg(PR->n, PR->m, PR->p, PR->C);
	y_best = y;
	f_best = funkcija_cilja(y, PR->n, PR->m, PR->p, PR->C);
	vremena(0) = merivreme(t1, sada());
	resenja(0) = f_best;
	pit = 0;
}
inline void insta::iteriraj() {
	vreme t1 = sada();
	y = konstrukcija_resenja_pohlepnim_alg(PR->n, PR->m, PR->p, PR->C);
	f = funkcija_cilja(y, PR->n, PR->m, PR->p, PR->C);
	ly = local_search(y, PR->n, PR->p, f, 0, 1, PR->C);
	lf = funkcija_cilja(ly, PR->n, PR->m, PR->p, PR->C);
	vremena(iter) = vremena(iter - 1) + merivreme(t1, sada());
	resenja(iter) = lf;
	if (lf < f_best) {
		y_best = ly;
		f_best = lf;
		pit = iter;
	}
	++iter;
}

uint64_t wm[4];
inline void pwm() {
	memset(wm, 0, 32);
	wm[0] |= 1;
	auto pp = [&](int lol) {
		int loli = lol >> 6;
		int loll = lol & 63;
		wm[loli] |= 1ull << loll;
	};
	pp(' ');
	pp('\t');
	pp('\n');
	pp('\r');
}
inline bool jeli(int lol) {
	int loli = lol >> 6;
	int loll = lol & 63;
	int l = wm[loli] >> loll;
	return l & 1;
}
inline double qeval(const char* ee, int l) {
	int c = -1;
	for (int i = 0; i < l; ++i)
		if (ee[i] == '/') {
			c = i;
			break;
		}
	if (c >= 0) {
		double a, b;
		sscanf(ee, "%lf", &a);
		sscanf(ee + c + 1, "%lf", &b);
		return a / b;
	}
	else {
		double a;
		sscanf(ee, "%lf", &a);
		return a;
	}
}
inline int intr(const char* ee, int l) {
	int x = 0;
	for (int i = 0; i < l; ++i)
		x = (x * 10) + ee[i] - '0';
	return x;
}
char patija[4096];
#ifdef __linux__
const char palim = '/';
#else
const char palim = '\\';
#endif
void* openauto(const char* filename, int64_t* size, int64_t limit = 0) {
	FILE* f = fopen(filename, "rb");
	if (f == nullptr) { *size = -1; return nullptr; }
	int64_t len = 0;
	try { len = std::filesystem::file_size(filename); }
	catch (...) { fclose(f); *size = -1; return nullptr; }
	if ((limit > 0) & (len > limit))len = limit;
	int64_t preso = 0;
	char* buc = (char*)malloc(len);
	while (preso < len) {
		int64_t isc = fread(buc + preso, 1, len - preso, f);
		if (isc < len - preso)
			if ((ferror(f) != 0) | (feof(f) != 0))
			{
				fclose(f); free(buc); *size = -1; return nullptr;
			}
		preso += isc;
	}
	fclose(f);
	*size = len;
	return buc;
}
int64_t save(const char* filename, void* buffer, int64_t length = 0, uint32_t blocksize = 4096) {
	FILE* f = fopen(filename, "wb");
	if (f == nullptr)return 0;
	if (length == 0)length = strlen((char*)buffer);
	int64_t tow = length / blocksize;
	int64_t r = fwrite(buffer, blocksize, tow, f);
	if (r < tow) { fclose(f); return r * blocksize; };
	int64_t ost = length % blocksize;
	if (ost) {
		r = fwrite(((char*)buffer) + (length - ost), ost, 1, f);
		fclose(f);
		return (length - ost) + r * ost;
	}
	fclose(f);
	return length;
}
struct heder {
	const char* name;
	int nlen;
	int goal;
	double q;
	inline heder(const char* _name, int _nlen, int _goal, double _q) :
		name(_name), nlen(_nlen), goal(_goal), q(_q) {}
};
// [folder sa primerima] [fajl sa listom primera] [default koeficijent] [output fajl]
const char* defa[5] = {
	nullptr,
	"D:\\Fakultet\\moje instance",
	"spisak.txt",
	"1.25",
	"auto_izvestaji.txt"
};
int main(int args, char** argv) {
	printf("de si\n");
#define RETURN(e) { printf("\nstisni enter da ugasis\n"); char xx; scanf("%c",&xx); return e; }
	const char** AA = (const char**)argv;
	if (args < 5) {
		// printf("invalid input\n");
		// RETURN(1)
		AA = defa;
	}

	arma_rng::set_seed_random();
	pwm();

	int pas = strlen(AA[1]);
	memcpy(patija, AA[1], pas);
	patija[pas++] = palim;

	int64_t ll;
	char* lista = (char*)openauto(AA[2], &ll);
	if (ll < 0) {
		printf("error\n");
		RETURN(2)
	}

	double stdq = qeval(AA[3], strlen(AA[3]));
	std::vector<heder> puri;

	int lul = 0;
	while (true) {
		/*
		   filename goal q |
		   filename goal |
		   |
		*/
		auto vadime = [&](int* poco, int* oco)->int {
			// cekaj crnog
			while (true) {
				if (lul == ll)return 2;
				if (lista[lul] == '\n') {
					++lul;
					return 1;
				}
				if (!jeli(lista[lul]))break;
				++lul;
			}
			*poco = lul;
			// cekaj belog
			while (true) {
				if (lul == ll) {
					*oco = lul;
					return 0;
				}
				if (jeli(lista[lul])) {
					*oco = lul;
					return 0;
				}
				++lul;
			}
		};
		/*
			0 - procito,
			1 - nije procito, novi red
			2 - nije procito, kraj fajla
		*/
		int poco_ime = 0;
		int oco_ime = 0;
		int poco_gol = 0;
		int oco_gol = 0;
		int poco_kf = -1;
		int oco_kf = -1;
		//citanje
		switch (vadime(&poco_ime, &oco_ime)) {
		case 0:	break;
		case 1: goto sledeci;
		case 2: goto napolje;
		}
		switch (vadime(&poco_gol, &oco_gol)) {
		case 0:	break;
		case 1: puri.push_back(heder(lista + poco_ime, oco_ime - poco_ime, -1, 0.0)); goto sledeci;
		case 2: puri.push_back(heder(lista + poco_ime, oco_ime - poco_ime, -1, 0.0)); goto napolje;
		}
		vadime(&poco_kf, &oco_kf);
		puri.push_back(heder(
			lista + poco_ime,
			oco_ime - poco_ime,
			intr(lista + poco_gol, oco_gol - poco_gol),
			poco_kf >= 0 ? qeval(lista + poco_kf, oco_kf - poco_kf) : stdq
		));
	sledeci:;
	}napolje:;

	int mfl = puri.size() * (450 + 50 * primer::INS);
	for (uint32_t i = 0; i < puri.size(); ++i)
		mfl += puri[i].nlen;
	char* fzc = (char*)malloc(mfl);
	if (fzc == nullptr) {
		printf("error\n");
		RETURN(3)
	}
	int sh = 0;

	primer P;
	for (uint32_t i = 0; i < puri.size(); ++i) {
		if (puri[i].goal < 0) {
			sh += sprintf(fzc + sh, "instance:   %2d, %.*s\n\n   bad input\n------------------------------------------------------------------------------------\n",
				i + 1, puri[i].nlen, puri[i].name);
			printf("%2d: bad input\n", i + 1);
			continue;
		}
		memcpy(patija + pas, puri[i].name, puri[i].nlen);
		memcpy(patija + pas + puri[i].nlen, ".txt", 5);
		if (P.ucitaj(patija, puri[i].q)) {
			sh += sprintf(fzc + sh, "instance:   %2d, %.*s\n\n   error\n------------------------------------------------------------------------------------\n",
				i + 1, puri[i].nlen, puri[i].name);
			printf("%2d: error\n", i + 1);
			continue;
		}
		printf("%2d: ", i + 1);
		fflush(stdout);
		P.resi(puri[i].goal);
		sh += P.print(i + 1, puri[i].name, puri[i].nlen, fzc + sh);
	}

	int osh = save(AA[4], fzc, sh);
	free(lista);
	if (osh == sh) {
		printf("sve je ispisano\n");
		free(fzc);
		RETURN(0)
	}

	int RE = 5 - (osh < 0);
	if (osh < 0) printf("nece da sacuva :D\n");
	else printf("ispisano je samo %d od %d bajtova\n", osh, sh);
	printf("\nza svaki slucaj, ovo je analiza:\n\n%.*s\n\n", sh, fzc);
	free(fzc);
	RETURN(RE)
}

