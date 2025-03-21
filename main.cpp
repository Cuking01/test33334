#include "CuSDK/basis.h"
#include "CuSDK/simd.h"

u3 pow_mod(u3 a,u3 b,u3 mod)
{
	u3 ans=1;
	while(b)
	{
		if(b&1)ans=(u4)ans*a%mod;
		a=(u4)a*a%mod;
		b>>=1;
	}
	return ans;
}

struct Mogo_F
{
	u3 mod;
	u3 modp;

	Mogo_F(u3 mod):mod(mod)
	{
		modp=calc_modp();
	}

	u3 calc_modp() const
	{
		return pow_mod(mod,(1ull<<51)-1,1ull<<52);
	}

	u3 calc_modp_newton() const
	{
		u3 x=1;
		for(int i=0;i<6;i++)
			x=x*2-x*x*mod;
		return x&((1ull<<52)-1);
	}


	void test_modp()
	{
		u3 x=calc_modp();
		printf("%llu %llu\n",x,x*mod%(1ull<<52));
	}

	void test_modp_newton()
	{
		u3 x=calc_modp_newton();
		printf("%llu %llu\n",x,x*mod%(1ull<<52));
	}

};

template<u2 n>
ALWAYS_INLINE void jmod(Pack_Ref<VU64x8,n> a,Pack_Ref<VU64x8,n> mod)
{
	Pack<VU64x8,n> tmp=a-mod;
	a=min(a,tmp);
}

Mogo_F_SIMD
{
	static inline const VU64x8 init_value=set(1ull,5ull,7ull,11ull,13ull,17ull,19ull,23ull);

	Pack<VU64x8,4> mod,modp;
	

	Mogo_F_SIMD(u3 mod_offset)
	{
		VU64x8 step=set1(24ull);
		VU64x8 offset=set1(mod_offset);
		mod[0]=init_value+offset;
		mod[1,2,3]=mod[0,1,2]+step;
		calc_modp();
	}

	void next()
	{
		VU64x8 step=set1(96ull);
		mod=mod+step;
		calc_modp();
	}

	void calc_modp()
	{
		VU64x8 zero=setzero(),two=set1(2ull);
		modp[0]=set1(1ull);
		modp[1,2,3]=modp[0];

		Pack<VU64x8,4> t;
		
		for(u2 i=0;i<6;i++)
		{
			t=madd52lo(two,modp,mod);
			modp=madd52lo(zero,modp,t);
		}
		modp=zero-modp;
	}

	template<bool need_jmod=true>
	void mul_mod(Pack_Ref<VU64x8,4> a,Pack_CRef<VU64x8,4> b)
	{
		Pack<VU64x8,4> xl,xh;
		VU64x8 zero=setzero();

		xl=madd52lo(zero,a,b);
		xh=madd52hi(zero,a,b);
		a=madd52lo(zero,xl,modp);
		a=madd52hi(xh,a,mod);

		if constexpr(need_jmod)
			jmod<4>(a,mod);
	}
};

int main()
{
	for(int i=1;i<=101;i+=2)
	{
		Mogo_F mf(i+(1ull<<50));

		mf.test_modp();
		mf.test_modp_newton();
	}
	

}
