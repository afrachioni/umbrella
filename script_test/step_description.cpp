#include <vector> 
#include <string>
#include "step_description.h"

StepDescription::StepDescription(){};

StepDescription::~StepDescription(){};

std::vector<std::string>* StepDescription::get_take_step_block() {
	return &take_step_block;
}

std::vector<std::string>* StepDescription::get_if_accept_block() {
	return &if_accept_block;
}

std::vector<std::string>* StepDescription::get_if_reject_block() {
	return &if_reject_block;
}

