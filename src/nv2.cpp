#include <stdio.h>
#include <cstdint>
#include <climits>
#include <cmath>
#include <thread>
#include <random>
#include <vector>
#include <homerus/annri.hpp>
double eval(const char* cmd)
{
	FILE* fp=fopen(".eval.cpp","w");
	fprintf(fp,"#include <stdio.h>\n");
	fprintf(fp,"#include <cmath>\n");
	fprintf(fp,"int main(void){\n");
	fprintf(fp,"using namespace std;\n");
	fprintf(fp,"\tprintf(\"%%lf\",(double)(%s));\n",cmd);
	fprintf(fp,"\treturn 0;\n");
	fprintf(fp,"}\n");
	fclose(fp);
	const char* mk="g++ -std=c++11 .eval.cpp -o .ev 2> /dev/null 1> /dev/null ; ./.ev > .ev.o";
	double val{(double)system(mk)};
	fp=fopen(".ev.o","r");
	if(fscanf(fp,"%lf",&val)!=1) val=0;
	fclose(fp);
	return val;
}
template <typename T>
struct ch_t
{
	T data[8];
	ch_t(){}
	T& operator[](size_t i){ return data[i]; }
};
struct cal_t
{
	double a,b;
	cal_t():a{1.},b{0.}{}
	cal_t(double a_,double b_):a{a_},b{b_}{}
	void set(double a_,double b_){ a=a_;b=b_;}
};
struct lin_bin_t
{
using index_t=size_t;
using value_t=double;
	index_t wid;
	value_t min,max;
	value_t ibin_;
	template <typename T1,typename T2,typename T3>
	lin_bin_t(T1 w,T2 n,T3 x)
		:wid{(index_t)(w)},min{(value_t)(n)},max{(value_t)(x)},
		ibin_{(value_t)(wid/(max-min))}
	{}
	index_t operator[](value_t val)
	{
		return (val<min||val>max)?index_t(-1):(index_t)((val-min)*ibin_);
	}
	value_t operator()(index_t n)
	{
		return min+n/ibin_;
	}
	value_t ibin(index_t ind)
	{
		return ibin_;
	}
};
struct large_bin_t
{
using index_t=size_t;
using value_t=double;
using small_bin_t=lin_bin_t;
	index_t      ste;
	index_t     mwid;
	value_t      min;
	value_t      max;
	value_t  logmin_;
	value_t ilogste_;
	index_t      wid;
	void default_(void)
	{
		ste=500;
		mwid=50;
		min  =1;
		max=5e7;
		logmin_=std::log(min);
		ilogste_=1./std::log(ste);
		wid=this->operator[](max);
	}
	template <typename T1,typename T2,typename T3,typename T4>
	large_bin_t(T1 st,T2 mw,T3 mi,T4 ma)
		:ste{(index_t)st},mwid{(index_t)mw},
		 min{(value_t)mi},max{(value_t)ma},logmin_{std::log(min)},ilogste_{1./std::log(value_t(ste))}
	{
		wid=this->operator[](max);
	}
	index_t operator[](value_t val)
	{
		uint32_t L=(uint32_t)((std::log(val)-logmin_)*ilogste_);
		value_t mi=min*std::pow(ste,L);
		value_t iw=mwid/((ste-1)*mi);
		return (val<min||val>max)?index_t(-1):(index_t)(L*mwid+(val-mi)*iw);
	}
	value_t operator()(index_t ind)
	{
		uint32_t N=uint32_t(ind/mwid);
		value_t p=min*std::pow(ste,N);
		return value_t(p+(ind-N*mwid)*p*(ste-1)/mwid);
	}
	value_t ibin(index_t ind)
	{
		uint32_t N=uint32_t(ind/mwid);
		value_t p=std::pow(ste,N);
		return mwid/(min*p*(ste-1));
	}
};
int main(int argc,char* argv[])
{
using namespace folklore::homerus;
using namespace std;
using ana::get_v;
using cal_ch_t=ch_t<cal_t>;
using bd_ch_t=ch_t< std::vector<size_t> >;
	char shm[]="/tmp/shm";
	char dmp[]=".map";
	char dcl[]=".cal";
	char dnv[]=".nv2";
	char stf[]="";
	char* rdf_fnm=shm;
	char* map_fnm=dmp;
	char* nv2_fnm=dnv;
	char* cal_fnm=dcl;
	char* set_fnm=stf;
	uint32_t eg[2]={0,8192+8192};
	uint32_t tg[2]={0,UINT_MAX-1};
	uint32_t ag[2]={0,UINT_MAX-1};
	auto mode = ana::offline;
	if(argc>1)
	{
		char onl[]="--online";
		char stg[]="--setting=";
		auto chk_=[](const char* s)->bool
		{
			for(;;)
			{
				if(*s=='\0') return false;
				if(*s==':' ) return true;
				++s;
			}
		};
		auto split_=[](char* s,uint32_t* buf)
		{
			if(sscanf(s+2,"%u:%u",buf,buf+1)!=2)
				printf("invalid syntax in \"%s\"\n",s);
			if(buf[0]>buf[1])
				printf("gate minimum seems to be smaller than gate maximum.\n");
		};
		for(int i=1;i<argc;++i)
		{
			char* ar=argv[i];
			if(chk_(ar))
			{
				switch(ar[0])
				{
					case 'e':
						split_(ar,eg);
						break;
					case 't':
						split_(ar,tg);
						break;
					case 'a':
						split_(ar,ag);
						break;
				}
			}
			else if(ar[0]!='-')
				rdf_fnm=ar;
			else if(ar[1]=='-'&&strncmp(ar,onl,strlen(onl))==0)
				mode=ana::online;
			else if(ar[1]=='-'&&strncmp(ar,stg,strlen(stg))==0)
			{
				size_t lstg=strlen(stg);
				size_t len=strlen(ar)-lstg;
				if(len==0) continue;
				char* fnm=(char*)malloc(len+1);
				strncpy(fnm,ar+lstg,len);
				set_fnm=fnm;
			}
			else
				printf("invalid syntax in \"%s\"\n",ar);
		}
	}
	auto replace_ch=[](char* b,char c,char r){
		for(size_t i=0;;++i)
		{
			if(b[i]==c) b[i]=r;
			if(b[i]=='\0') break;
		}
	};
	{
		FILE *fp=fopen(set_fnm,"r");
		if(fp)
		{
			char buf[512];
			while(fgets(buf,512,fp))
			{
				replace_ch(buf,'\n','\0');
				replace_ch(buf,'=',' ');
				char cmd[512];
				char arg[512];
				char* tmp;
				tmp=(char*)malloc(strlen(arg));
				strcpy(tmp,arg);
				if(strcmp(cmd,"nv2")==0) nv2_fnm=tmp;
				if(strcmp(cmd,"cal")==0) cal_fnm=tmp;
				if(strcmp(cmd,"map")==0) map_fnm=tmp;
				if(strcmp(cmd,"rdf")==0) rdf_fnm=tmp;
			}
		}
	}
	homerus<ana::annri> hom(rdf_fnm,mode);
	if(!hom)
	{
		printf("cannot open file: %s\n",rdf_fnm);
		return 0;
	}
	double e_min{0};
	double e_max{8192*2};
	double t_ste{10.};
	double t_min{1};
	double t_max{8192*2};
	uint32_t e_wid{8192*2};
	uint32_t t_wid{8192*2};
	{
		FILE* fp=fopen(nv2_fnm,"r");
		if(fp)
		{
			auto count_ch=[](char* str,char c)->uint32_t{
				uint32_t v{0};
				for(size_t i=0;str[i]!='\n';++i) if(str[i]==c) ++v;
				return v;
			};
			char buf[512];
			while(fgets(buf,512,fp)!=nullptr)
			{
				if(buf[0]=='#') continue;
				if(strncmp(buf,"e_bin",strlen("e_bin"))==0)
				{
					if(count_ch(buf,',')!=2) continue;
					char* b[3]{buf+6,buf+6,buf+6};
					for(size_t i=6,j=1;buf[i]!='\n';++i) if(buf[i]==','){buf[i]='\0'; b[j++]=buf+i+1;}
					e_wid=(uint32_t)eval(b[0]);
					e_min=eval(b[1]);
					e_max=eval(b[2]);
				}
				if(strncmp(buf,"t_bin",strlen("t_bin"))==0)
				{
					if(count_ch(buf,',')!=3) continue;
					char* b[4]{buf+6,buf+6,buf+6,buf+6};
					for(size_t i=6,j=1;buf[i]!='\n';++i) if(buf[i]==','){buf[i]='\0'; b[j++]=buf+i+1;}
					t_ste=(uint32_t)eval(b[0]);
					t_wid=(uint32_t)eval(b[1]);
					t_min=eval(b[2]);
					t_max=eval(b[3]);
				}
			}
			fclose(fp);
		}
	}
	lin_bin_t   ebin{e_wid,e_min,e_max};
	large_bin_t tbin{t_ste,t_wid,t_min,t_max};
	vector<cal_ch_t> ecal;
	vector<cal_ch_t> tcal;
	{
		FILE *fp=fopen(cal_fnm,"r");
		bool is_e{false};
		bool is_t{false};
		if(fp)
		{
			char bf[512];
			size_t ln{0};
			while(fgets(bf,512,fp))
			{
				if(bf[0]=='#'||bf[0]=='\n') { ++ln; continue; }
				for(size_t i=0;;++i) if(bf[i]=='#'){ bf[i]='\n'; break; }else if(bf[i]=='\n') break;
				char* tmp=bf;
				for(size_t i=0;;++i) if(bf[i]!=' '&&bf[i]!='\t'){ tmp=bf+i; break; }
				if(tmp[0]=='#'||tmp[0]=='\n') { ++ln; continue; }
				if(tmp[0]=='[')
				{
					if(strncmp(tmp,"[[ecal]]",strlen("[[ecal]]"))==0){ is_e=true; is_t=false; }
					if(strncmp(tmp,"[[tcal]]",strlen("[[tcal]]"))==0){ is_t=true; is_e=false; }
				}
				else
				{
					size_t dt[2]; double val[2];
					if(sscanf(tmp,"%zu %zu %lf %lf",dt,dt+1,val,val+1)!=4) 
					{
						printf("syntax error at line #%zu in file %s\n",ln++,cal_fnm);
						continue;
					}
					if(!is_e && !is_t) continue;
					vector<cal_ch_t>* v=(is_e?&ecal:&tcal);
					if((dt[0]+1)>v->size()) for(size_t i=v->size();i<dt[0]+1;++i) v->push_back(cal_ch_t{});
					if(dt[1]<8)((*v)[dt[0]])[dt[1]].set(val[0],val[1]);
				}
				++ln;
			}
		}
	}
	uint32_t max_id{0};
	vector<bd_ch_t> map;
	{
		FILE *fp=fopen(map_fnm,"r");
		if(fp)
		{
			char bf[512];
			size_t ln{0};
			while(fgets(bf,512,fp))
			{
				if(bf[0]=='#'||bf[0]=='\n') { ++ln; continue; }
				char a[512];
				uint32_t id{0},bd{0};
				if(sscanf(bf,"%u %u %s",&id,&bd,a)!=3)
				{
					printf("syntax error at line #%zu in file %s\n",ln++,map_fnm);
					continue;
				}
				if(bd+1>map.size())
					for(size_t i=map.size();i<bd+1;++i) map.push_back(bd_ch_t{});
				max_id=id>max_id?id:max_id;
				if(strcmp(a,"*")==0)
				{
					for(size_t i=0;i<8;++i)
						map[bd][i].push_back(id);
					continue;
				}
				auto ch=atoi(a);
				if(ch>=0&&ch<8) map[bd][ch].push_back(id);
			}
		}
	}
	if(max_id==0) max_id=31;
	++max_id;
	auto print_info=[&](FILE* fp)
	{
		fprintf(fp,"# file = %s\n",rdf_fnm);
		fprintf(fp,"# energy gate = %u - %u\n",eg[0],eg[1]);
		fprintf(fp,"# time gate   = %u - %u\n",tg[0],tg[1]);
		fprintf(fp,"# ebin        = %zu,%lf,%lf\n",ebin.wid,ebin.min,ebin.max);
		fprintf(fp,"# tbin        = %zu,%zu,%lf,%lf\n",tbin.ste,tbin.mwid,tbin.min,tbin.max);
	};
	print_info(stdout);
	bool start{false};
	bool update{false};
	char ctr{'a'};
	random_device  dev;
	mt19937 eng(dev());
	uniform_real_distribution<> dis(-0.5,0.5);
	double* eh=(double*)malloc(sizeof(double)*ebin.wid*(2+max_id));
	double* th=(double*)malloc(sizeof(double)*tbin.wid*(2+max_id));
	memset(eh,0,sizeof(double)*ebin.wid*(2+max_id));
	memset(th,0,sizeof(double)*tbin.wid*(2+max_id));
	for(size_t i=0;i<ebin.wid;++i)
	{
		eh[i]         =ebin(i);
		eh[ebin.wid+i]=1./ebin.ibin(i);
	}
	for(size_t i=0;i<tbin.wid;++i)
	{
		th[i]         =tbin(i);
		th[tbin.wid+i]=1./tbin.ibin(i);
	}
	auto update_ehs=[&](const char* ehn,bool fl)
	{
		FILE *fp=fopen(ehn,"w");
		if(fl) print_info(fp);
		for(uint32_t i=0;i<ebin.wid;++i)
		{
			fprintf(fp,"%lf\t%lf\t",eh[i],eh[ebin.wid+i]);
			double in=1./eh[ebin.wid+i];
			for(uint32_t j=2;j<2+max_id;++j)
				fprintf(fp,"%g\t",eh[ebin.wid*j+i]*in);
			fprintf(fp,"\n");
		}
		fclose(fp);
	};
	auto update_ths=[&](const char* thn,bool fl)
	{
		FILE *fp=fopen(thn,"w");
		if(fl) print_info(fp);
		for(uint32_t i=0;i<tbin.wid;++i)
		{
			fprintf(fp,"%lf\t%lf\t",th[i],th[tbin.wid+i]);
			double in=1./th[tbin.wid+i];
			for(uint32_t j=2;j<2+max_id;++j)
				fprintf(fp,"%g\t",th[tbin.wid*j+i]*in);
			fprintf(fp,"\n");
		}
		fclose(fp);
	};
	auto fill=[&](vector<size_t> l,double e,double t){
		if(l.empty()) return;
		if(!((e>eg[0]&&e<eg[1])&&(t>tg[0]&&t<tg[1]))) return;
		if(e>ebin.min&&e<ebin.max)
			for(auto id:l) eh[ebin.wid*(2+id)+ebin[e]]+=1;
		if(t>tbin.min&&t<tbin.max)
			for(auto id:l) th[tbin.wid*(2+id)+tbin[t]]+=1;
	};
	auto update_eval=[&](const char* ehn,const char* thn,bool fl)
	{
		update_ehs(ehn,fl);
		update_ths(thn,fl);
		update=false;
	};
	std::thread update_([&]()
	{
		while(!start) sleep(1);
		for(;;)
		{
			if(update) update_eval("ehs.tmp","ths.tmp",false);
			sleep(2);
		}
	});
	update_.detach();
	std::thread key([&]()
	{
		while(!start) sleep(1);
		printf("press \'q<Enter>\' to quit program\n");
		for(;;)
		{
			if(scanf("%c",&ctr)!=1) continue;
			usleep(1000);
		}
	});
	key.detach();
	size_t cnt{0};
	for(auto data:hom)
	{
		start=true;
		++cnt;
		if(cnt<ag[0]) continue;
		if(cnt>ag[1]) break;
		auto e = get_v<ana::energy >(data);
		auto t = get_v<ana::time   >(data);
		auto b = get_v<ana::board  >(data); 
		auto c = get_v<ana::channel>(data); 
		if(map.size()<=b) continue;
		fill(
			map[b][c],
			ecal.empty()?e:ecal[b][c].a*(dis(eng)+e)+ecal[b][c].b,
			tcal.empty()?t:tcal[b][c].a*(dis(eng)+t)+tcal[b][c].b
		);
		if(ctr=='q') break;
		update=true;
	}
	update_eval("ehs.dat","ths.dat",true);
	return 0;
}
