#include <iostream>
#include <thread>
#include <mutex>
#include <vector>
#include <atomic>
using namespace std;
std::function<void()> func;
mutex mt;
std::condition_variable cv;
struct test
{
	test(int a=1, int b=1) :a(a), b(b) {}
	int a, b;
	void fun()
	{
		func = [&] {
			unique_lock<mutex> loc(mt);
			cv.wait(loc);
			cout << a << ' ' << b << endl;
		};
	}
};
vector<thread>pool;
int thread_count = 8;
int now_count = 0;
int main()
{
	for (int i = 1; i <= 7; ++i)
	{
		pool.push_back(thread([&] (int out){
			{
				std::unique_lock<mutex> loc(mt);
				thread_count--;
				if (!thread_count)
					now_count++,cv.notify_all();
				else
					now_count++,cv.wait(loc);
			}
			cout << now_count << endl;
				//while (1)
				//	cout << out << endl;
		},i));
	}
	{
		std::unique_lock<mutex> loc(mt);
		thread_count--;
		if (!thread_count)
			now_count++, cv.notify_all();
		else
			now_count++, cv.wait(loc);
	}
	cout << now_count << endl;
	for (int i = 0; i < pool.size(); ++i)
		pool[i].join();
	//while (1)
	//	cout << 8 << endl;
	//test now;
	//now.fun();
	//std::thread t(func);
	//cv.notify_all();
	////std::thread t(&test::func,&now);
	//t.join();
    std::cout << "Hello World!\n"<<endl;
}

