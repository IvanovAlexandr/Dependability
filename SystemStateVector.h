/*
  @file SystemStateVector.h
  @brief Implementation computing of reliability
  @author Alexandr Ivanov (ivanov.alexandr.1995@gmail.com)
*/

#ifndef DEPENDABILITY_SYSTEM_STATE_VECTOR_H
#define DEPENDABILITY_SYSTEM_STATE_VECTOR_H


#include <cstdint>
#include <vector>

class SystemStateVector {
public:

    unsigned long getSizeVcc1zero();
    unsigned long getSizeVcc2zero();
    unsigned long getSizeVcc3zero();
    unsigned long getSizeVcc4zero();

    void generateVcc();
    void generateVirtualSystemStateVector();

    std::vector<double> getPassiveFaultTolerance();
    std::vector<double> getActiveFaultTolerance();


private:

    void createVirtualSystemStateVector(std::vector<std::pair<uint32_t, uint32_t>> &vcc);


    std::vector<std::pair<uint32_t, uint32_t>> m_vcc1zero;
    std::vector<std::pair<uint32_t, uint32_t>> m_vcc2zero;
    std::vector<std::pair<uint32_t, uint32_t>> m_vcc3zero;
    std::vector<std::pair<uint32_t, uint32_t>> m_vcc4zero;

};


#endif //DEPENDABILITY_SYSTEM_STATE_VECTOR_H
