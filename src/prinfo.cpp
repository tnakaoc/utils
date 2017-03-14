#include <stdio.h>
#include <string.h>
#include <stdlib.h>
void remove_cr(char* ar,size_t max)
{
	for(size_t i=0;i<max;++i)
	{
		if(ar[i]=='\0') break;
		if(ar[i]=='\n')
		{
			ar[i]='\0';
			break;
		}
	}
}
int main(int argc,char* argv[])
{
	if(argc==1)
	{
		printf("usage %s [filename]\n",argv[0]);
		return 0;
	}
	FILE *fp=fopen(argv[1],"rb");
	if(!fp)
	{
		printf("file open error (%s)\n",argv[1]);
		return 0;
	}
	char buffer[512];
	if(fread(buffer,512,1,fp)==1)
	{
		auto pinfo=[&](){
			char tmp[256];
			size_t offset[3]={64,128,256};
			for(auto off:offset)
			{
				strncpy(tmp,buffer+off,off);
				remove_cr(tmp,off);
				printf("%s\n",tmp);
			}
			printf("\n");
		};
		printf("header info:\n");
		pinfo();
		fseek(fp,0,SEEK_END);
		fseek(fp,-512,SEEK_CUR);
		fread(buffer,512,1,fp);
		printf("ender info:\n");
		pinfo();
	}
	else
	{
		printf("invalid filetype\n");
	}
	fclose(fp);
	return 0;
}

