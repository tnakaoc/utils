#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <iostream>
#include <cmath>
#include <include/homerus/annri.hpp>
#include <include/large_bin.hpp>
template <uint32_t S>
struct hist
{
	double sum{0};
	double cen{0};
	double sig{0};
	uint32_t cnt{0};
	const size_t G0=3500000;
	const size_t G1=8800000;
	hist()
	{
		clear();
	}
	uint32_t buf[S];
	void clear()
	{
		memset(buf,0,sizeof(buf));
	}
	void insert(size_t s,size_t g)
	{
		if(s<S && g>G0 && g<G1) ++buf[s];
	}
	uint32_t operator[](size_t s) const
	{
		return buf[s];
	}
	constexpr uint32_t size() const
	{
		return S;
	}
	void stat(size_t beg,size_t en)
	{
		double xsum{0};
		sum = 0;
		cen=0;
		sig=0;
		cnt=0;
		if(beg>=S) return;
		if(en<beg || en>S) en=S;
		double isum{0};
		for(size_t i=beg;i<en;++i)
		{
			isum+=i;
			sum +=buf[i];
			xsum+=buf[i]*i;
		}
		cen=xsum/sum;
		xsum=0;
		for(size_t i=beg;i<en;++i)
		{
			double R=buf[i]*(cen-i)*(cen-i);
			xsum+=R;
		}
		sig=std::sqrt(xsum)/sum;
	}
};
int main(int argc,char* argv[])
{
using namespace folklore::homerus;
using ana::get_v;
using ana::energy;
	if(argc==1 || std::string(argv[1])=="-h" || std::string(argv[1])=="--help")
	{
		printf("\n");
		printf("%s : extended viewer with peak shift correction.(spesialized for Li glass analysys.)\n",argv[0]);
		printf("usage: %s [filename] ([OPTIONS])\n",argv[0]);
		printf("\n");
		printf("OPTIONS:\n");
		printf("\teg   : energy gate\n");
		printf("\t\tf.e. eg10:1000:20:1200 -> gate 10 to 1000 for GS20, 20 to 1200 for GS30\n");
		printf("\ttg   : timing gate\n");
		printf("\t\tsame as nv or viewer (common values for GS20 and GS30)\n");
		printf("\tedge : energy gate edge to find initial peak positions.\n");
		printf("\t\tf.e. edge300:1000:200:1200 \n\t\t\t-> search peak in ch300 to ch1000 for GS20, ch200 to ch1200 for GS30\n");
		printf("\tpc   : peak center\n");
		printf("\n");
		printf("command sample:\n");
		printf("\t%s file.rdf eg2500:4000:2000:4000 tg10000:2e7 edge3000:5000:2000:6000\n",argv[0]);
		printf("\tmeans,\n");
	   	printf("\t\tanalyse \"file.rdf\" with energy gate 2500-4000 for GS20, 2000-4000 for GS30\n");
		printf("\t\tpeak searching window is 3000-5000 for GS20, 2000-6000 for GS30\n");
		printf("\n");
		return 0;
	}
	uint32_t edge[2][2]={{2000,5000},{2000,4000}};
	uint32_t eg[2][2]={{0,8192*2},{0,8192*2}};
	uint32_t tg[2]={0,UINT32_MAX-1};
	uint32_t ag[2]={0,0};
	uint32_t shot[2]={15,15};
	double center{3500};
	homerus<ana::annri> hom(argv[1]);
	if(!hom)
	{
		printf("file open error (%s)\n",argv[1]);
		return 0;
	}
	if(argc!=2)
	{
		double tmpv[4]{0};
		for(int i=2;i<argc;++i)
		{
			switch(argv[i][0])
			{
				case 'e':
					if(argv[i][1]=='g')
					{
						sscanf(argv[i],"eg%lf:%lf:%lf:%lf",tmpv,tmpv+1,tmpv+2,tmpv+3);
						for(uint32_t k=0;k<2;++k) for(uint32_t l=0;l<2;++l)
							if(tmpv[2*k]>0 &&tmpv[2*k]<tmpv[2*k+1] && tmpv[2*k+1]<8192*2) eg[k][l]=(uint32_t)tmpv[2*k+l];
					}
					else if(argv[i][1]=='d')
					{
						sscanf(argv[i],"edge%lf:%lf:%lf:%lf",tmpv,tmpv+1,tmpv+2,tmpv+3);
						for(uint32_t k=0;k<2;++k) for(uint32_t l=0;l<2;++l)
							if(tmpv[2*k]>0 &&tmpv[2*k]<tmpv[2*k+1] && tmpv[2*k+1]<8192*2) edge[k][l]=(uint32_t)tmpv[2*k+l];
					}
					break;
				case 't':
					sscanf(argv[i],"tg%lf:%lf",tmpv,tmpv+1);
					if(tmpv[0]>0 && tmpv[0]<tmpv[1] && tmpv[1]<UINT32_MAX-1)
					{
						tg[0]=(uint32_t)tmpv[0];tg[1]=(uint32_t)tmpv[1];
					}
					break;
				case 'a':
					sscanf(argv[i],"ag%lf:%lf",tmpv,tmpv+1);
					if(tmpv[0]>=0 && tmpv[0]<tmpv[1])
					{
						ag[0]=(uint32_t)tmpv[0];
						ag[1]=(uint32_t)tmpv[1];
					}
					break;
				case 'p':
					sscanf(argv[i],"pc%lf",tmpv);
					if(tmpv[0]>0 && tmpv[0]<5000) center=tmpv[0];
					break;
				default:
					break;
			}
		}
	}
	uint32_t count[2]={1,1};
	folklore::large_bin<500,2> bin;
	hist<8192*2> h[2];
	const uint32_t HSIZE=9;
	const size_t FOSIZE=8192;
	const size_t FOBIN =8192*8;
	uint64_t *ehist[HSIZE];
	uint64_t *thist[HSIZE];
	uint32_t fobuf[FOSIZE];
	uint32_t prevt[2]{0,0};
	const uint32_t E_max{8192*2};
	const uint32_t T_max=bin[UINT32_MAX-1];
	printf("# T bin max = %u\n",T_max);
	for(size_t i=0;i<HSIZE;++i)
	{
		ehist[i]=(uint64_t*)malloc(sizeof(uint64_t)*E_max);
		thist[i]=(uint64_t*)malloc(sizeof(uint64_t)*T_max);
		memset(ehist[i],0,sizeof(uint64_t)*E_max);
		memset(thist[i],0,sizeof(uint64_t)*T_max);
	}
	memset(fobuf,0,sizeof(fobuf));
	std::vector<uint32_t> evs[2];
	std::vector<uint32_t> tvs[2];
	for(uint32_t i=0;i<2;++i)
	{
		evs[i].clear();
		tvs[i].clear();
		evs[i].push_back(0);
		tvs[i].push_back(0);
	}
	FILE* fp=fopen("/dev/urandom","rb");
	auto get_rand=[&](void)->double
	{
		char c;
		fread(&c,sizeof(char),1,fp);
		return ((double)c)/256;
	};
	uint32_t ch2id[8]={0,2,1,3,3,3,3,3};
	uint32_t ac{0};
	for(auto data:hom)
	{
		auto bd=get_v<ana::board>(data);
		if(bd!=0) continue;
		auto cc=get_v<ana::channel>(data);
		if(cc>=8) continue;
		auto id=ch2id[cc];
		if(id==3) continue;
		auto e=get_v<energy>(data);
		if(id==2)
		{
			e=e>E_max?1:e;
			++ehist[id][e];
			++thist[id][1];
			++ac;
			if(ag[0]>=ac) continue;
			if(ag[1]!=0 && ag[1]<=ac) break;
			continue;
		}
		auto tm=get_v<ana::time>(data);
		if(tm<prevt[id])
		{
			++count[id];
			if(count[id]>=shot[id])
			{
				count[id]=0;
				h[id].stat(edge[id][0],edge[id][1]);
				if(false)
				if(h[id].cen>edge[id][0] && h[id].cen<edge[id][1])
				{
					uint32_t from = (uint32_t)(h[id].cen-h[id].sig*2);
					uint32_t to   = (uint32_t)(h[id].cen+h[id].sig*2);
					h[id].stat(from,to);
				}
				auto rat=center/h[id].cen;
				size_t previd=0;
				auto end=evs[id].size();
				for(size_t j=1;j<end;++j)
				{
					auto   ev=evs[id][j];
					auto   tv=tvs[id][j];
					double ec=ev*rat+get_rand();
					uint32_t iec=(uint32_t)ec;
					if(iec>0&&iec<E_max) ++ehist[id][iec];
					if(ev>0 && ev<E_max) ++ehist[id+5][ev];
					auto bn=bin[tv];
					if(bn>0 && bn<T_max)
					{
						++thist[id][bn];
						thist[id+7][bn]=tv;
					}
					if((tv<tvs[id][j-1])||j==end-1)
					{
						if((tvs[id][j-1]>12000000)||j==end-1)
						{
							uint32_t lc{0};
							uint32_t beg=previd;
							for(auto II=beg;II<j-1;++II)
							{
								auto TT=tvs[id][II];
								if(TT<10000000) ++lc;
								else break;
								beg=II;
							}
							beg=lc>100?beg+1:previd;
							for(auto II=beg;II<j-1;++II)
							{
								auto TT=tvs[id][II];
								auto EE=evs[id][II];
								double ecc=EE*rat+get_rand();
								if(ecc>eg[id][0]&&ecc<eg[id][1]&&TT>tg[0]&&TT<tg[1])
								{
									auto Bn=bin[TT];
									++thist[id+5][Bn];
									uint32_t bnn{(uint32_t(TT/FOBIN))};
									if(id==0&&bnn<FOSIZE)
									{
										++fobuf[bnn];
									}
								}
							}
						}
						previd=j;
					}
					if(ec>eg[id][0]&&ec<eg[id][1]&&tv>tg[0]&&tv<tg[1])
					{
						if(iec>0&&iec<E_max) ++ehist[id+3][iec];
						if(bn>0 && bn<T_max) ++thist[id+3][bn];
					}
		 		}
		 		h[id].clear();
		 		evs[id].clear();
		 		tvs[id].clear();
				evs[id].push_back(0);
				tvs[id].push_back(0);
				prevt[id]=0;
			}
		}
		evs[id].push_back(e);
		tvs[id].push_back(tm);
		h[id].insert(e,tm);
		prevt[id]=tm;
	}
	const char *efname="ehs.env";
	const char *tfname="ths.env";
	printf("# scanf finished for %s\n",argv[1]);
	printf("# output files are %s and %s\n",efname,tfname);
	auto print_fileheader=[&](FILE* fp){
		fprintf(fp,"# file: %s\n",argv[1]);
		fprintf(fp,"# analyse           : %u-%u events\n",ag[0],ag[1]);
		fprintf(fp,"# energy gate(GS20) : %u-%u ch\n",eg[0][0],eg[0][1]);
		fprintf(fp,"# energy gate(GS30) : %u-%u ch\n",eg[1][0],eg[1][1]);
		fprintf(fp,"# time gate         : %u-%u ch\n",tg[0],tg[1]);
		fprintf(fp,"# peak center       : %lf ch\n",center);
		fprintf(fp,"# window(GS20)      : %u-%u ch\n",edge[0][0],edge[0][1]);
		fprintf(fp,"# window(GS30)      : %u-%u ch\n",edge[1][0],edge[1][1]);
	};
	{
		FILE* fp=fopen(efname,"w");
		print_fileheader(fp);
		fprintf(fp,"#ch GS20 GS30 Trigger GS20(gated) GS30(gated) GS20(raw) GS30(raw)\n");
		for(size_t i=0;i<E_max;++i)
		{
			fprintf(fp,"%zu\t",i);
			for(size_t j=0;j<HSIZE;++j)
				fprintf(fp,"%lu\t",ehist[j][i]);
			fprintf(fp,"\n");
		}
		fclose(fp);
	}
	{
		FILE* fp=fopen("fo.env","w");
		print_fileheader(fp);
		fprintf(fp,"#ch GS20\n");
		for(size_t j=0;j<FOSIZE;++j)
			fprintf(fp,"%zu\t%lf\n",j*FOBIN,(double)fobuf[j]/FOBIN);
		fclose(fp);
	}
	{
		FILE* fp=fopen(tfname,"w");
		print_fileheader(fp);
		fprintf(fp,"#ch GS20 GS30 Trigger GS20(gated) GS30(gated) GS20(frame) GS30(frame) GS20(bootstrap T) GS30(bootstrap T) BinWidth\n");
		for(size_t i=0;i<T_max;++i)
		{
			auto tval = thist[7][i];
			uint32_t wid  = bin(tval);
			if(wid==0) continue;
			double inv[HSIZE]={1./wid,1./wid,1./wid,1./wid,1./wid,1./wid,1./wid,1,1};
			auto tch = bin.round(tval);
			fprintf(fp,"%u\t",(uint32_t)tch);
			for(size_t j=0;j<HSIZE;++j)
				fprintf(fp,"%lf\t",(double)(thist[j][i])*inv[j]);
			size_t tpos=(size_t)(tval/FOBIN);
			fprintf(fp,"%lf\t",tpos<FOSIZE?((double)(fobuf[tpos]))/FOBIN:0.);
			fprintf(fp,"%u\n",wid);
		}
		fclose(fp);
	}
	return 0;
}
