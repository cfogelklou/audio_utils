#define USE_CURL 0

#include "src/roughtime_request.hpp"
#include "src/roughtime_parse.hpp"

#if (USE_CURL > 0)
#include <curl/curl.h>
#include <curl/easy.h>
#endif
#include <gtest/gtest.h>
#include <gmock/gmock.h>
using namespace testing;

#ifdef WIN32
#include <Windows.h>
#define usleep(x) Sleep(x/1000)
#endif

#if (USE_CURL > 0)
TEST(TestCurl, curl_1) {
  CURL * curl = curl_easy_init();
  EXPECT_NE(nullptr, curl);
  curl_easy_cleanup(curl);  
}
#endif


static void stupidRandom(uint8_t *buf, int cnt) {
  for (int i = 0; i < cnt; i++) {
    buf[i] = rand() % 255;
  }
}

TEST(TestRt, UnpaddedRequest){
  sstring req;
  uint8_t nonce[64];
  stupidRandom(nonce, sizeof(nonce));
  RoughTime::GenerateRequest(req, nonce, sizeof(nonce));
  ASSERT_EQ(req.length(), 80);
  const uint8_t *pr = req.data();
  EXPECT_EQ(pr[0], 2); // 2 tags
  EXPECT_EQ(pr[1], 0);
  EXPECT_EQ(pr[2], 0);
  EXPECT_EQ(pr[3], 0);
  EXPECT_EQ(pr[4], 64); // Second tag has offset 64
  EXPECT_EQ(pr[5], 0);
  EXPECT_EQ(pr[6], 0);
  EXPECT_EQ(pr[7], 0);
  EXPECT_EQ(pr[8], 'N'); // First data (at offset zero) has tag nonc
  EXPECT_EQ(pr[9], 'O');
  EXPECT_EQ(pr[10], 'N');
  EXPECT_EQ(pr[11], 'C');
  EXPECT_EQ(pr[12], 'P'); // Second data (at offset 64) has tag PAD
  EXPECT_EQ(pr[13], 'A');
  EXPECT_EQ(pr[14], 'D');
  EXPECT_EQ(pr[15], 0xff);

}

TEST(TestRt, PaddedRequest){
  sstring req;
  uint8_t nonce[64];
  stupidRandom(nonce, sizeof(nonce));
  RoughTime::GenerateRequest(req, nonce, sizeof(nonce));
  ASSERT_EQ(req.length(), 80);
  RoughTime::PadRequest(req, req);
  ASSERT_EQ(req.length(), 1024);
  const uint8_t *pr = req.data();
  EXPECT_EQ(pr[0], 2); // 2 tags
  EXPECT_EQ(pr[1], 0);
  EXPECT_EQ(pr[2], 0);
  EXPECT_EQ(pr[3], 0);
  EXPECT_EQ(pr[4], 64); // Second tag has offset 64
  EXPECT_EQ(pr[5], 0);
  EXPECT_EQ(pr[6], 0);
  EXPECT_EQ(pr[7], 0);
  EXPECT_EQ(pr[8], 'N'); // First data (at offset zero) has tag nonc
  EXPECT_EQ(pr[9], 'O');
  EXPECT_EQ(pr[10], 'N');
  EXPECT_EQ(pr[11], 'C');
  EXPECT_EQ(pr[12], 'P'); // Second data (at offset 64) has tag PAD
  EXPECT_EQ(pr[13], 'A');
  EXPECT_EQ(pr[14], 'D');
  EXPECT_EQ(pr[15], 0xff);

  EXPECT_EQ(pr[1023], 0xff);
}

TEST(TestRt, RoughtimeParse) {
  const uint8_t nonce[] = { 5,83,77,186,201,145,69,152,83,175,116,202,170,224,249,101,238,22,56,197,147,36,14,207,95,96,87,228,209,115,210,207,11,73,211,160,163,96,120,215,35,70,189,113,239,168,19,66,87,253,196,183,159,189,59,195,190,205,149,160,113,25,228,80 };

  EXPECT_EQ(sizeof(nonce) , 64);

  const uint8_t packet[] = { 2,0,0,0,64,0,0,0,78,79,78,67,80,65,68,255,5,83,77,186,201,145,69,152,83,175,116,202,170,224,249,101,238,22,56,197,147,36,14,207,95,96,87,228,209,115,210,207,11,73,211,160,163,96,120,215,35,70,189,113,239,168,19,66,87,253,196,183,159,189,59,195,190,205,149,160,113,25,228,80,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 };
  EXPECT_EQ(sizeof(packet) , 1024);

  const uint8_t answer[] = { 5,0,0,0,64,0,0,0,64,0,0,0,164,0,0,0,60,1,0,0,83,73,71,0,80,65,84,72,83,82,69,80,67,69,82,84,73,78,68,88,73,155,232,86,166,133,223,183,220,49,82,201,230,189,142,241,222,152,21,121,97,17,75,234,111,156,147,90,106,7,48,125,177,46,67,105,6,109,248,234,16,74,219,210,232,133,135,104,159,76,117,3,78,81,60,143,146,168,252,219,110,156,119,8,3,0,0,0,4,0,0,0,12,0,0,0,82,65,68,73,77,73,68,80,82,79,79,84,64,66,15,0,216,217,100,144,127,131,5,0,141,188,245,106,82,96,124,212,119,88,228,116,35,184,129,193,14,215,109,9,121,26,34,251,219,96,97,97,30,98,245,79,94,255,156,34,30,190,95,254,163,232,173,188,236,78,119,119,126,255,231,14,106,5,217,89,187,226,103,228,225,232,87,90,2,0,0,0,64,0,0,0,83,73,71,0,68,69,76,69,136,199,40,177,0,54,169,60,144,114,166,1,236,0,51,68,40,41,2,246,32,74,27,154,240,238,32,180,205,203,87,222,183,178,218,164,235,190,114,157,246,106,48,40,136,132,197,142,165,56,245,236,230,253,140,206,151,136,82,39,134,200,155,15,3,0,0,0,32,0,0,0,40,0,0,0,80,85,66,75,77,73,78,84,77,65,88,84,3,197,172,255,171,147,174,153,132,224,166,216,116,206,1,95,50,69,25,215,47,249,61,190,241,62,229,10,115,56,255,216,48,213,14,9,115,131,5,0,48,53,230,38,135,131,5,0,0,0,0,0 };
  EXPECT_EQ(sizeof(answer), 360);

  const uint8_t pubkey[] = { 128, 62, 183, 133, 40, 247, 73, 196, 190, 194, 227, 158, 26, 187, 155, 94, 90, 183, 228, 221, 92, 228, 182, 242, 253, 47, 147, 236, 195, 83, 143, 26 };
  EXPECT_EQ(sizeof(pubkey), 32);

  const uint64_t midpoint = 1551958790167000;
  const uint64_t radius = 1000000;

  RoughTime::ParseOutT times;
  EXPECT_EQ(midpoint, RoughTime::ParseToMicroseconds(pubkey, nonce, answer, sizeof(answer), &times));

  EXPECT_EQ(times.radius, radius);
  EXPECT_EQ(times.midpoint, midpoint);

}

TEST(TestRt, RoughtimeSend) {
  const uint8_t nonce[] = { 242,255,68,93,235,98,11,53,248,11,11,42,185,146,154,114,6,153,211,252,204,238,171,235,78,0,245,40,116,124,115,181,176,55,95,43,116,251,88,40,60,117,141,163,203,151,226,90,202,211,43,252,217,50,60,192,151,61,7,26,216,137,133,110 };
  const uint8_t answer[] = { 5, 0, 0, 0, 64, 0, 0, 0, 64, 0, 0, 0, 164, 0, 0, 0, 60, 1, 0, 0, 83, 73, 71, 0, 80, 65, 84, 72, 83, 82, 69, 80, 67, 69, 82, 84, 73, 78, 68, 88, 217, 138, 43, 155, 232, 101, 19, 17, 204, 180, 107, 70, 130, 81, 248, 45, 0, 180, 198, 96, 109, 59, 211, 118, 172, 90, 96, 15, 1, 227, 147, 25, 70, 246, 16, 182, 20, 114, 41, 95, 117, 116, 168, 97, 207, 147, 252, 125, 121, 181, 99, 11, 120, 197, 17, 135, 206, 129, 41, 236, 214, 154, 193, 2, 3, 0, 0, 0, 4, 0, 0, 0, 12, 0, 0, 0, 82, 65, 68, 73, 77, 73, 68, 80, 82, 79, 79, 84, 64, 66, 15, 0, 160, 150, 146, 91, 127, 131, 5, 0, 192, 126, 226, 147, 9, 106, 35, 252, 9, 234, 110, 238, 23, 13, 26, 171, 10, 104, 181, 70, 53, 66, 250, 121, 48, 157, 11, 207, 220, 143, 163, 125, 39, 148, 28, 148, 69, 7, 66, 230, 77, 124, 157, 112, 139, 110, 169, 100, 148, 55, 172, 27, 186, 235, 206, 102, 229, 3, 140, 186, 200, 11, 55, 247, 2, 0, 0, 0, 64, 0, 0, 0, 83, 73, 71, 0, 68, 69, 76, 69, 200, 136, 163, 241, 194, 76, 82, 232, 237, 194, 224, 130, 205, 120, 66, 192, 91, 170, 105, 135, 125, 29, 135, 90, 103, 186, 255, 61, 22, 226, 99, 57, 254, 167, 174, 15, 121, 36, 116, 8, 77, 64, 150, 167, 188, 57, 161, 32, 122, 207, 62, 151, 146, 5, 212, 177, 185, 89, 252, 131, 231, 191, 27, 7, 3, 0, 0, 0, 32, 0, 0, 0, 40, 0, 0, 0, 80, 85, 66, 75, 77, 73, 78, 84, 77, 65, 88, 84, 178, 44, 182, 6, 18, 160, 216, 188, 147, 171, 123, 88, 36, 97, 220, 213, 191, 42, 77, 203, 24, 205, 201, 211, 26, 179, 250, 213, 180, 205, 155, 172, 152, 37, 242, 169, 119, 131, 5, 0, 152, 133, 201, 199, 139, 131, 5, 0, 0, 0, 0, 0 };
  const uint8_t pubkey[] = { 128, 62, 183, 133, 40, 247, 73, 196, 190, 194, 227, 158, 26, 187, 155, 94, 90, 183, 228, 221, 92, 228, 182, 242, 253, 47, 147, 236, 195, 83, 143, 26 };
  const uint64_t midpoint = 1551957903972000;
  const uint64_t radius = 1000000;
}

#if 1
#include "mini_socket.hpp"

TEST(TestRt, SendRequest){
  sstring req;
  uint8_t nonce[64];
  stupidRandom(nonce, sizeof(nonce));
  RoughTime::GenerateRequest(req, nonce, sizeof(nonce));
  ASSERT_EQ(req.length(), 80);
  RoughTime::PadRequest(req, req);
  ASSERT_EQ(req.length(), 1024);
  const uint8_t *pr = req.data();
  EXPECT_EQ(pr[0], 2); // 2 tags
  EXPECT_EQ(pr[1], 0);
  EXPECT_EQ(pr[2], 0);
  EXPECT_EQ(pr[3], 0);
  EXPECT_EQ(pr[4], 64); // Second tag has offset 64
  EXPECT_EQ(pr[5], 0);
  EXPECT_EQ(pr[6], 0);
  EXPECT_EQ(pr[7], 0);
  EXPECT_EQ(pr[8], 'N'); // First data (at offset zero) has tag nonc
  EXPECT_EQ(pr[9], 'O');
  EXPECT_EQ(pr[10], 'N');
  EXPECT_EQ(pr[11], 'C');
  EXPECT_EQ(pr[12], 'P'); // Second data (at offset 64) has tag PAD
  EXPECT_EQ(pr[13], 'A');
  EXPECT_EQ(pr[14], 'D');
  EXPECT_EQ(pr[15], 0xff);
  EXPECT_EQ(pr[1023], 0xff);


  QueueBase &p = CreateUdpClient("roughtime.cloudflare.com", 2002);
  p.Write(req.data(), req.length());
  size_t bytes = 0;
  int tries = 0;
  while ((tries < 1000) && (0 == bytes)){
    usleep(100000);
    bytes = p.GetReadReady();
    tries++;
  }

  EXPECT_GT(bytes, 0u);
  if (bytes) {
    uint8_t buf[1024];
    auto amtRead = p.Read(buf, bytes);
    EXPECT_EQ(amtRead, bytes);

    const uint8_t pubkey[] = { 128, 62, 183, 133, 40, 247, 73, 196, 190, 194, 227, 158, 26, 187, 155, 94, 90, 183, 228, 221, 92, 228, 182, 242, 253, 47, 147, 236, 195, 83, 143, 26 };
    RoughTime::ParseOutT times;
    
    uint64_t midpoint = RoughTime::ParseToMicroseconds(pubkey, nonce, buf, amtRead, &times);
    midpoint /= (1000ull*1000ull);

    std::time_t now = std::time(0);
    uint64_t compareTo =  (uint64_t)now;
    auto diff = fabs(compareTo - midpoint);
    EXPECT_LE(diff,  2000);

  }
  
  DeleteUdpClient(&p);
     
}
#endif


extern "C" {
  // This is a terrible random number generator, but only for unit testing.
  void randombytes_buf(uint8_t * const buf, const uint32_t size) {
    for (uint32_t i = 0; i < size; i++) {
      buf[i] = rand() & 0xff;
    }
  }

  void sodium_misuse(void) {
    assert(0);
  }
}

int main(int argc, char** argv){
  
  // The following line must be executed to initialize Google Mock
  // (and Google Test) before running the tests.
  ::testing::InitGoogleMock(&argc, argv);
  const int gtest_rval = RUN_ALL_TESTS();
  
  return gtest_rval;
}