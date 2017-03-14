#ifndef FOLKLORE_LARGE_BIN_HPP
#define FOLKLORE_LARGE_BIN_HPP
namespace folklore
{
	namespace histgram
	{
		struct xaxis_scale_log{};
		struct xaxis_scale_linear{};
	}
	template <long Seg_,long Step_,typename Tr_=folklore::histgram::xaxis_scale_log>
	struct large_bin
	{
		typedef long value_t;
		struct bin_t { bin_t(value_t b,value_t w):bin(b),width(w){} value_t bin; value_t width; };
		bin_t ch2bin_t(unsigned long ch)
		{
			if(ch<Seg_) return bin_t(0,0);
			int digit=(int)(std::log(ch/Seg_)/std::log(Step_));
			int sbin =(int)std::pow(Step_,digit+1);
			return bin_t((Step_-1)*digit*Seg_/Step_+ (int)((ch-sbin*Seg_/Step_)/sbin) ,sbin);
		}
		value_t ch2bin(unsigned long ch)
		{
			return ch2bin_t(ch).bin;
		}
		value_t ch2width(unsigned long ch)
		{
			return ch2bin_t(ch).width;
		}
		value_t operator[](unsigned long ch)
		{
			return ch2bin_t(ch).bin;
		}
		value_t operator()(unsigned long ch)
		{
			return ch2bin_t(ch).width;
		}
		value_t round(unsigned long ch)
		{
			if(ch<Seg_) return 0;
			int digit=(int)(std::log(ch/Seg_)/std::log(Step_));
			int sbin =(int)std::pow(Step_,digit+1);
			value_t offset=Seg_*sbin/Step_;
			return offset+sbin*(int)((ch-offset)/sbin);
		}
	};
}
#endif // FOLKLORE_LARGE_BIN_HPP
