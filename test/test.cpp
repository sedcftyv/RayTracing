#include <iostream>
#include <thread>
using namespace std;
std::function<void()> func;
struct test
{
	test(int a=1, int b=1) :a(a), b(b) {}
	int a, b;
	void fun()
	{
		func = [&] {
			cout << a << ' ' << b << endl;
		};
	}
};
int main()
{
	test now;
	now.fun();
	std::thread t(func);
	//std::thread t(&test::func,&now);
	t.join();
    std::cout << "Hello World!\n"<<endl;
}

