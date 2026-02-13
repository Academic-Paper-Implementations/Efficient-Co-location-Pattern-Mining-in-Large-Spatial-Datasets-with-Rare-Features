/**
 * @file constants.h
 * @brief Named constants for the joinless colocation mining algorithm
 * 
 * This file defines all magic numbers as named constants to improve
 * code readability and maintainability.
 */

#pragma once

namespace Constants {
    // Epsilon values for numerical stability
    constexpr double EPSILON_SMALL = 1e-9;      ///< Small epsilon for division by zero protection
    constexpr double EPSILON_DELTA = 1e-9;      ///< Delta threshold for rare intensity calculations
    
    // Progress reporting
    constexpr double MAX_PROGRESS_PERCENT = 95.0;  ///< Maximum progress percentage before completion
    
    // Default values
    constexpr double DEFAULT_MIN_PREVALENCE = 0.6;  ///< Default minimum prevalence threshold
    constexpr double DEFAULT_MIN_COND_PROB = 0.5;   ///< Default minimum conditional probability
    constexpr double DEFAULT_NEIGHBOR_DISTANCE = 5.0; ///< Default neighbor distance threshold
}
