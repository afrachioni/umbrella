#include <vector> 
#include <string>
#include "step_description.h"

StepDescription::StepDescription(){};

StepDescription::~StepDescription(){};

std::vector<std::string>* StepDescription::get_take_step_block() {
	return &take_step_block;
}

