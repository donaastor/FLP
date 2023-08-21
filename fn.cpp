#include <cstdio>
#include <cstdlib>
#include <cstdint>
#include <climits>
#include <cstring>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <chrono>
#include <filesystem>
#include <vector>
#ifndef __linux__
#pragma warning(disable:4996)
#endif

/*
	matrica C: n(.i.) puta m(.j.)
		C[i,j] je udaljenost i-tog coveka od j-tog mesta
		cuva se tako da je j monotono, a i nije
	resenje je odabir p mesta (predstavljeno binarnim vektorim duzine m)
		kada svakog coveka spojimo sa njemu najblizem mestu,
		racunamo razliku najvise i najmanje spojenog mesta
*/

enum{ // strategy
	FIRST_IMPROV=0,
	BEST_IMPROV=1
};
const int STRATEGY = FIRST_IMPROV
;

// funkcijice

uint64_t randin=1;
inline void randround(){
	randin^=randin>>12;
	randin^=randin<<25;
	randin^=randin>>27;
	randin*=0x2545f4914f6cdd1dull;
}
inline int vmax(int *v,int m){
	int r=INT_MIN;
	for(int j=0;j<m;++j)
		if(v[j]>r)
			r=v[j];
	return r;
}
inline int vmin(int *v,int m){
	int r=INT_MAX;
	for(int j=0;j<m;++j)
		if(v[j]<r)
			r=v[j];
	return r;
}
inline int randint(int m){ // nasumicno bira ceo broj iz intervala [0,m)
	randround();
	return (randin>>24)%m;
}
inline int randcoef(int *v,int m){ // nasumicno bira jednu koordinatu vektora
	return v[randint(m)];
}
int funkcija_cilja(char *res,int n,int m,int p,int *C){
	int najblm[n];
	int tabl[n];
	int koll[m];
	int j=0;
	while((j<m)&(res[j]^1))++j; // nadji prvog keca
	for(int i=0;i<n;i++)najblm[i]=j;
	memcpy(tabl,C+j*n,4*n); // iskopiraj udaljenosti
	for(;++j<m;){ // idi po ostalim kecevima
		if(res[j]==0)continue;
		int *tl=C+j*n;
		for(int i=0;i<n;++i) // izminuj udaljenosti
			if(tl[i]<tabl[i]){
				tabl[i]=tl[i];
				najblm[i]=j;
			}
	}
	memset(koll,0,4*m);
	for(int i=0;i<n;++i) // izbroj povezanosti
		koll[najblm[i]]++;
	int smin=INT_MAX,smax=INT_MIN;
	for(j=0;j<m;++j) // nadji najmanju i najvecu
		if(res[j]){
			if(koll[j]<smin) smin=koll[j];
			if(koll[j]>smax) smax=koll[j];
		}
	return smax-smin;
}
inline void tumbani(int n,int m,int *C,int *Cl){ // Cl je duzine m
	memset(Cl,0,4*m);
	for(int j=0,k=0;j<m;++j)
		for(int i=0;i<n;++i,++k)
			if(C[k]!=j+1)++Cl[j];
}
void konstrukcija_resenja_pohlepnim_alg(char *res,int n,int m,int p,int *C){
	memset(res,0,m);
	res[randint(m)]=1;
	int Cl[m]; tumbani(n,m,C,Cl);
	int Rcl[m];
	for(int i=1;i<p;++i){
		int a=INT_MAX,b=INT_MIN;
		for(int j=0;j<m;++j){
			if(res[j])continue;
			if(Cl[j]<a) a=Cl[j];
			if(Cl[j]>b) b=Cl[j];
		}
		int l=0,c=ceil(a+0.5*(b-a));
		for(int j=0;j<m;++j)
			if((Cl[j]>=c)&(res[j]^1))
				Rcl[l++]=j;
		res[randcoef(Rcl,l)]=1;
	}
}
inline void inversion(char *res,int poc,int kraj){
	for(int i=poc,j=kraj;i<j;++i,--j){
		char t=res[i];
		res[i]=res[j];
		res[j]=t;
	}
}
template<const int strategy>
int local_search(char *res,int f,int n,int m,int p,int *C){
	char dres[m];
	int najx=0,najy=0,najf=f;
	dalje:
	for(int y=1;y<m;++y)for(int x=0;x<y;++x){
		memcpy(dres,res,m);
		inversion(dres,x,y);
		int nf=funkcija_cilja(dres,n,m,p,C);
		if(nf<najf){
			najf=nf;
			if constexpr(strategy==FIRST_IMPROV){
				memcpy(res,dres,m);
				goto dalje;
			}
			else{
				najx=x;
				najy=y;
			}
		}
	}
	if constexpr(strategy==FIRST_IMPROV)return f;
	else if(najf==f)return f;
	inversion(res,najx,najy);
	f=najf;
	goto dalje;
}

// merenje vremena

using namespace std::chrono;
typedef steady_clock::time_point vreme;
inline vreme sada(){
	return steady_clock::now();
}
inline double merivreme(vreme poc,vreme kraj){
	return duration_cast<duration<double>>(kraj-poc).count();
}

// resavanje

struct primer;
struct insta{
	#define PRnmpC PR->n,PR->m,PR->p,PR->C
	primer *PR;
	char *res;
	int f,it,best_it;
	double v_best;
	double v_end;
	inline void pisivreme(const char *txt);
	inline void pisivreme(const char *txt,int ix);
	inline void pocni(primer *group);
	inline void skoncaj();
	inline void iteriraj();
};
struct primer{
	int n,m,p;
	int maxiter,it;
	int *C;
	constexpr static int INS = 15
;	insta ie[INS];
	inline int ucitaj(const char *fname,double mik){
		FILE *file=fopen(fname,"r");
		if(file==nullptr)
			return -1;
		fscanf(file,"%d",&m);
		fscanf(file,"%d",&n);
		fscanf(file,"%d",&p);
		maxiter=n*mik;
		C=(int*)malloc(4*n*m);
		for(int i=0;i<n;++i)for(int j=0;j<m;++j)
			fscanf(file,"%d",C+j*n+i);
		fclose(file);
		return 0;
	}
	inline void smakni(){
		free(C);
	}
	inline int maxres(){
		int najg=ie[0].f;
		for(int i=1;i<INS;++i)
			if(ie[i].f>najg)
				najg=ie[i].f;
		return najg;
	}
	inline int minres(){
		int najb=ie[0].f;
		for(int i=1;i<INS;++i)
			if(ie[i].f<najb)
				najb=ie[i].f;
		return najb;
	}
	inline void resi(int goal){
		int ispiso=0;
		static const char *tems="                                                                                ";
		static const char *zevs="||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||";
		for(int i=0;i<INS;++i)ie[i].pocni(this);
		int res=maxres();
		if(res<=goal){
			it=0;
			for(int i=0;i<INS;++i)ie[i].skoncaj();
			printf("%s done\n",tems);
			return;
		}
		bool ng=false;
		int zavi=maxiter;
		it=-1;
		for(int j=1;j<=zavi;++j){
			for(int i=0;i<INS;++i)
				ie[i].iteriraj();
			int nisp=(80*j)/maxiter;
			printf("%.*s",nisp-ispiso,zevs);
			fflush(stdout);
			ispiso=nisp;
			res=maxres();
			if((ng==false)&(res<=goal)){
				ng=true;
				zavi=j*2;
				if(zavi>maxiter)zavi=maxiter;
				it=zavi;
			}
		}
		for(int i=0;i<INS;++i)ie[i].skoncaj();
		printf(" done\n");
	}
	inline double ukvreme(){
		double x=0;
		for(int i=0;i<INS;++i)
			x+=ie[i].v_end;
		return x;
	}
	inline double najvreme(){
		double x=0;
		for(int i=0;i<INS;++i)
			x+=ie[i].v_best;
		return x;
	}
	inline int print(int ord,const char *name,int nlen,char *os){
		int sh=0;
		sh+=sprintf(os+sh,"instance:   %2d, %.*s\nm,n,p:      %d %d %d\n",ord,nlen,name,m,n,p);
		sh+=sprintf(os+sh,"\n       sol_i         t_i      ttot_i       gap_i    best_sol           t        ttot        agap       sigma        iter\n");
		int best_sol=minres();
		double avgt=najvreme(),avgttot=ukvreme(),ssol=0,qsol=0;
		for(int i=0;i<INS;++i)ssol+=ie[i].f;
		for(int i=0;i<INS;++i)qsol+=ie[i].f*ie[i].f;
		avgt/=INS;avgttot/=INS;ssol/=INS;qsol/=INS;
		const double kf=100.0/best_sol;
		qsol=kf*sqrt(qsol-ssol*ssol);
		ssol=kf*(ssol-best_sol);
		char gaps[13];gaps[12]=0;
		auto cgap=[&](int ind)->const char*{
			double g=kf*(ie[ind].f-best_sol);
			if(g!=0) sprintf(gaps,"%12.2f",g);
			else memcpy(gaps,"           -",12);
			return gaps;
		};
		sh+=sprintf(os+sh,"%12d%12.5f%12.5f%s%12d%12.5f%12.5f%12.5f%12.5f%12d\n",
			ie[0].f,ie[0].v_best,ie[0].v_end,cgap(0),
			best_sol,avgt,avgttot,ssol,qsol,it);
		for(int i=1;i<INS;++i)
			sh+=sprintf(os+sh,"%12d%12.5f%12.5f%s\n",
				ie[i].f,ie[i].v_best,ie[i].v_end,cgap(i));
		sh+=sprintf(os+sh,"------------------------------------------------------------------------------------------------------------------------\n");
		return sh;
	}
};
inline void insta::pocni(primer *group){
	PR=group;
	res=(char*)malloc(PR->m);
	vreme v_poc=sada();
	konstrukcija_resenja_pohlepnim_alg(res,PRnmpC);
	f=funkcija_cilja(res,PRnmpC);
	v_best=v_end=merivreme(v_poc,sada());
	it=best_it=1;
}
inline void insta::skoncaj(){
	free(res);
}
inline void insta::iteriraj(){
	vreme v_poc=sada();
	char tres[PR->m];
	konstrukcija_resenja_pohlepnim_alg(tres,PRnmpC);
	int nf=local_search<STRATEGY>(tres,funkcija_cilja(tres,PRnmpC),PRnmpC);
	v_end+=merivreme(v_poc,sada());
	++it;
	if(nf<f){
		memcpy(res,tres,PR->m);
		f=nf;
		best_it=it;
		v_best=v_end;
	}
}

// mucenje sa fajlovima

uint64_t wm[4];
inline void pwm(){
	memset(wm,0,32);
	wm[0]|=1;
	auto pp=[&](int lol){
		int loli=lol>>6;
		int loll=lol&63;
		wm[loli]|=1ull<<loll;
	};
	pp(' ');
	pp('\t');
	pp('\n');
	pp('\r');
}
inline bool jeli(int lol){
	int loli=lol>>6;
	int loll=lol&63;
	int l=wm[loli]>>loll;
	return l&1;
}
inline double qeval(const char *ee,int l){
	int c=-1;
	for(int i=0;i<l;++i)
		if(ee[i]=='/'){
			c=i;
			break;
		}
	if(c>=0){
		double a,b;
		sscanf(ee,"%lf",&a);
		sscanf(ee+c+1,"%lf",&b);
		return a/b;
	}
	else{
		double a;
		sscanf(ee,"%lf",&a);
		return a;
	}
}
inline int intr(const char *ee,int l){
	int x=0;
	for(int i=0;i<l;++i)
		x=(x*10)+ee[i]-'0';
	return x;
}
char patija[4096];
#ifdef __linux__
const char palim='/';
#else
const char palim='\\';
#endif
void *openauto(const char *filename,int64_t *size,int64_t limit=0){
	FILE *f=fopen(filename,"rb");
	if(f==nullptr){*size=-1;return nullptr;}
	int64_t len=0;
	try{len=std::filesystem::file_size(filename);}
	catch(...){fclose(f);*size=-1;return nullptr;}
	if((limit>0)&(len>limit))len=limit;
	int64_t preso=0;
	char *buc=(char*)malloc(len);
	while(preso<len){
		int64_t isc=fread(buc+preso,1,len-preso,f);
		if(isc<len-preso)
			if((ferror(f)!=0)|(feof(f)!=0))
				{fclose(f);free(buc);*size=-1;return nullptr;}
		preso+=isc;
	}
	fclose(f);
	*size=len;
	return buc;
}
int64_t save(const char *filename,void *buffer,int64_t length=0,uint32_t blocksize=4096){
	FILE *f=fopen(filename,"wb");
	if(f==nullptr)return 0;
	if(length==0)length=strlen((char*)buffer);
	int64_t tow=length/blocksize;
	int64_t r=fwrite(buffer,blocksize,tow,f);
	if(r<tow){fclose(f);return r*blocksize;};
	int64_t ost=length%blocksize;
	if(ost){
		r=fwrite(((char*)buffer)+(length-ost),ost,1,f);
		fclose(f);
		return (length-ost)+r*ost;
	}
	fclose(f);
	return length;
}
struct heder{
	const char *name;
	int nlen;
	int goal;
	double q;
	inline heder(const char *_name,int _nlen,int _goal,double _q):
		name(_name),nlen(_nlen),goal(_goal),q(_q){}
};
/// [folder sa primerima] [fajl sa listom primera] [default koeficijent] [output fajl]
const char *defa[5]={
	nullptr,
	"D:\\Fakultet\\moje instance",
	"spisak.txt",
	"1.25",
	"auto_izvestaji.txt"
};
int main(int args,char **argv){
	#define RETURN(e) { printf("\nstisni enter da ugasis\n"); getchar(); return e; }
	const char **AA=(const char**)argv;
	if(args<5){
		printf("invalid input\n");
		// RETURN(1)
		AA=defa;
	}
	
	randin=time(0);
	randin+=randin==0;
	pwm();
	
	int pas=strlen(AA[1]);
	memcpy(patija,AA[1],pas);
	patija[pas++]=palim;
	
	int64_t ll;
	char *lista=(char*)openauto(AA[2],&ll);
	if(ll<0){
		printf("error\n");
		RETURN(2)
	}
	
	double stdq=qeval(AA[3],strlen(AA[3]));
	std::vector<heder> puri;
	
	int lul=0;
	while(true){
		/*
		   filename goal q |
		   filename goal |
		   |
		*/
		auto vadime=[&](int *poco,int *oco)->int{
			// cekaj crnog
			while(true){
				if(lul==ll)return 2;
				if(lista[lul]=='\n'){
					++lul;
					return 1;
				}
				if(!jeli(lista[lul]))break;
				++lul;
			}
			*poco=lul;
			// cekaj belog
			while(true){
				if(lul==ll){
					*oco=lul;
					return 0;
				}
				if(jeli(lista[lul])){
					*oco=lul;
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
		int poco_ime=0;
		int oco_ime=0;
		int poco_gol=0;
		int oco_gol=0;
		int poco_kf=-1;
		int oco_kf=-1;
		//citanje
		switch(vadime(&poco_ime,&oco_ime)){
			case 0:	break;
			case 1: goto sledeci;
			case 2: goto napolje;
		}
		switch(vadime(&poco_gol,&oco_gol)){
			case 0:	break;
			case 1: puri.push_back(heder(lista+poco_ime,oco_ime-poco_ime,-1,0.0));goto sledeci;
			case 2: puri.push_back(heder(lista+poco_ime,oco_ime-poco_ime,-1,0.0));goto napolje;
		}
		vadime(&poco_kf,&oco_kf);
		puri.push_back(heder(
			lista+poco_ime,
			oco_ime-poco_ime,
			intr(lista+poco_gol,oco_gol-poco_gol),
			poco_kf>=0?qeval(lista+poco_kf,oco_kf-poco_kf):stdq
		));
	sledeci:
	}napolje:
	
	int mfl=puri.size()*(450+50*primer::INS);
	for(uint32_t i=0;i<puri.size();++i)
		mfl+=puri[i].nlen;
	char *fzc=(char*)malloc(mfl);
	if(fzc==nullptr){
		printf("error\n");
		RETURN(3)
	}
	int sh=0;
	
	primer P;
	for(uint32_t i=0;i<puri.size();++i){
		if(puri[i].goal<0){
			sh+=sprintf(fzc+sh,"instance:   %2d, %.*s\n\n   bad input\n------------------------------------------------------------------------------------\n",
				i+1,puri[i].nlen,puri[i].name);
			printf("%2d: bad input\n",i+1);
			continue;
		}
		memcpy(patija+pas,puri[i].name,puri[i].nlen);
		memcpy(patija+pas+puri[i].nlen,".txt",5);
		if(P.ucitaj(patija,puri[i].q)){
			sh+=sprintf(fzc+sh,"instance:   %2d, %.*s\n\n   error\n------------------------------------------------------------------------------------\n",
				i+1,puri[i].nlen,puri[i].name);
			printf("%2d: error\n",i+1);
			continue;
		}
		printf("%2d: ",i+1);
		fflush(stdout);
		P.resi(puri[i].goal);
		P.smakni();
		sh+=P.print(i+1,puri[i].name,puri[i].nlen,fzc+sh);
	}
	
	int osh=save(AA[4],fzc,sh);
	free(lista);
	if(osh==sh){
		printf("sve je ispisano\n");
		free(fzc);
		RETURN(0)
	}
	
	int RE=5-(osh<0);
	if(osh<0) printf("nece da sacuva :D\n");
	else printf("ispisano je samo %d od %d bajtova\n",osh,sh);
	printf("\nza svaki slucaj, ovo je analiza:\n\n%.*s\n\n",sh,fzc);
	free(fzc);
	RETURN(RE)
}

