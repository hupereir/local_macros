#include <random>

void test_matching()
{
  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_int_distribution<uint64_t> distrib(0, 0xFFFFFFFFFFU);

  for( int i=0; i < 1000; ++i )
  {
    const uint64_t first = distrib(gen)&0xFFFFFFFFFFU;
    const uint64_t second = distrib(gen)&0xFFFFFU;

    const uint64_t diff_1 = ((first&0xFFFFFU) - second)&0xFFFFFU;
    const uint64_t diff_2 = (first- second)&0xFFFFFU;

    std::cout << "0x" << std::hex << diff_1 << " 0x" << diff_2 << std::dec << std::endl;
  }
}
