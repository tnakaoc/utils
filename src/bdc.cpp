#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <set>
#include <homerus/annri.hpp>
int main(int argc,char* argv[])
{
using namespace folklore::homerus;
using namespace std;
using ana::get_v;
	if(argc==1)
	{
		printf("usage: %s [filename]\n",argv[0]);
		return 0;
	}
	set<uint32_t> bdlist;
	homerus<ana::annri> hom(argv[1]);
	for(auto data:hom)
	{
		bdlist.insert(get_v<ana::board>(data));
	}
	for(auto i:bdlist)
	{
		printf("%u\n",i);
	}
	return 0;
}

