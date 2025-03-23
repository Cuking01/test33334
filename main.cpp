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

struct Tester
{
	static inline const VU64x8 init_value=set(1ull,7ull,11ull,13ull,17ull,19ull,23ull,29ull);

	Pack<VU64x8,4> mod,modp;  //模数和蒙哥马利模乘参数
	Pack<VU64x8,4> x;
	VU64x8 flag,ans;

	u3 mod_offset;
	u3 i;
	u3 num;

	ALWAYS_INLINE Tester(u3 mod_offset,u3 num):mod_offset(mod_offset),i(0),num(num)
	{
		VU64x8 step=set1(30ull);
		VU64x8 offset=set1(mod_offset);
		mod[0]=init_value+offset;
		mod[1,2,3]=mod[0,1,2]+step;
	}

	ALWAYS_INLINE void pow333()
	{
		Pack<VU64x8,4> t;
		for(int i=0;i<27;i++)
		{
			t=x;
			mul_mod<false>(t,t);
			mul_mod<false>(x,t);
		}
	}

	template<bool check_algorithm=false>
	ALWAYS_INLINE void check()
	{
		calc_modp();
		enter_mogo();
		x[2].print("enter");
		pow333();
		x[2].print("pow333");
		leave_mogo();
		x[2].print("leave");

		VU64x8 four=set1(4ull),zero=set1(0ull),mask=set1((1ull<<52)-1);
		x=x+four;
		jmod<4>(x,mod);
		x=x&mask;

		if constexpr(check_algorithm)
		{
			static constexpr u3 offset[8]={1ull,7ull,11ull,13ull,17ull,19ull,23ull,29ull};

			auto calc_r=[](u3 mod)
			{
				u3 ans=pow_mod(3,7625597484987,mod);
				ans=(ans+4)%mod;
				return ans;
			};

			alignas(64) u3 tmp[32];
			x.store(tmp);
			for(int j=0;j<32;j++)
				if(tmp[j]!=calc_r(mod_offset+i+offset[j%8]+(j/8)*30))
				{
					puts("Error.");
					printf("%llu %llu\n",mod_offset,i);
					x[0].print("x[0]");
					x[1].print("x[1]");
					x[2].print("x[2]");
					x[3].print("x[3]");

					for(int j=0;j<32;j++)
					{
						printf("%llu ",calc_r(mod_offset+i+offset[j%8]+(j/8)*30));
						if(j%8==7)puts("");
					}
					exit(0);
				}
		}

		x[0,1]=min(x[0,1],x[2,3]);
		flag=min(x[0],x[1]);

	}

	//测试素性
	template<bool check_algorithm=false>
	ALWAYS_INLINE bool test()
	{
		VU64x8 step=set1(120ull);
		while(i<num)
		{
			check<check_algorithm>();
			mod=mod+step;
			i+=120;
		}

		alignas(64) u3 tmp[8];
		flag.store(tmp);
		for(int i=0;i<8;i++)
			if(tmp[i]==0)
				return false;
		return true;
	}

	ALWAYS_INLINE void calc_modp()
	{
		VU64x8 zero,two=set1(2ull);
		zero.setzero();
		modp[0]=set1(1ull);
		modp[1,2,3]=modp[0];

		Pack<VU64x8,4> t;
		
		for(u2 i=0;i<6;i++)
		{
			t=madd52lo(two,modp,mod);
			modp=madd52lo(zero,modp,t);
		}
	}

	ALWAYS_INLINE void enter_mogo()
	{
		u3 k=(1ull<<52)/(mod_offset+i+119ull);
		VU64x8 vk=set1(k);
		VU64x8 zero;
		zero.setzero();

		modp[2].print("modp");

		x=madd52lo(zero,vk,mod);
		x=zero-x;   //求出1的蒙哥马利域的值，但可能会差若干个p
		x[2].print("mogo_1'\n");
		mul_mod(x,x); //求1的平方，依然是1，但是值的范围会压缩到[0,p)内
		x[2].print("mogo_1\n");
		Pack<VU64x8,4> y;
		y=x+x;
		y=y+x; //求出3在蒙哥马利域的值。
		mul_mod(x,y); //1*3，得到3并把范围压缩到[0,p);
		x[2].print("mogo_3\n");
	}

	ALWAYS_INLINE void leave_mogo()
	{
		VU64x8 zero;
		zero.setzero();
		VU64x8 one=set1(1ull);

		x=madd52lo(zero,x,modp);
		x=madd52hi(one,x,mod);

		jmod<4>(x,mod);
	}

	template<bool need_jmod=true>
	ALWAYS_INLINE void mul_mod(Pack_Ref<VU64x8,4> a,Pack_CRef<VU64x8,4> b)
	{
		Pack<VU64x8,4> xl,xh;
		VU64x8 zero,one=set1(1ull);
		zero.setzero();

		xl=madd52lo(zero,a,b);
		xh=madd52hi(one,a,b);
		a=madd52lo(zero,xl,modp);
		a=madd52hi(xh,a,mod);

		if constexpr(need_jmod)
			jmod<4>(a,mod);
	}
};



#include <thread>
#include <condition_variable>

bool test_for(u3 l,u3 block_size,u3 block_num,u3 thread_num)
{
	std::mutex mtx;
	std::condition_variable cv;
	std::atomic<u3> cnt=0;
	std::atomic<u3> block_id=0;
	std::atomic<bool> passed=true;
	std::vector<std::thread> work_threads;
	work_threads.reserve(thread_num);

	auto worker=[&mtx,&cv,&cnt,&block_id,&passed,l,block_size,block_num]()
	{
		while(1)
		{
			u3 cur_id=block_id++;
			if(cur_id>=block_num)break;

			Tester tester(l+cur_id*block_size,block_size);
			if(!tester.test())
			{
				passed=false;
				break;
			}
		}

		cnt++;
		cv.notify_all();
	};

	for(u3 i=0;i<thread_num;i++)
		work_threads.emplace_back(worker);
	
	std::unique_lock lock(mtx);
	cv.wait(lock,[&]{return cnt==thread_num||passed==false;});

	for(auto&thread:work_threads)
		thread.join();

	return passed;
}

void check_algorithm(u3 l)
{
	// for(u3 i=0;i<100000000000000ull;i+=100000000000ull)
	// {
	// 	Tester tester(l+i,120*100);
	// 	tester.test<true>();
	// }

	Tester tester(52301200000000ull+9240,120);
	tester.test<true>();
}

int main()
{
	u3 l=120*10000000;

	check_algorithm(l);

	int st=clock();

	bool ret=test_for(l,120ull*833334,1000,32);

	int ed=clock();

	printf("time=%d\n",ed-st);

	if(ret)
	{
		puts("may be a prime.");
	}
	else
	{
		puts("not a prime.");
	}
}
