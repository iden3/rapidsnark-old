#include <gtest/gtest.h>
#include <alt_bn128.hpp>

AltBn128::Engine E;
 
int main(int argc, char** argv) { 
    ::testing::InitGoogleTest(&argc, argv);
    
    return RUN_ALL_TESTS();
}
