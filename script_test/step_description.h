

class StepDescription {
	public:
		StepDescription();
		~StepDescription();
		std::vector<std::string>* get_take_step_block();
		std::vector<std::string>* get_if_accept_block();
		std::vector<std::string>* get_if_reject_block();

		void take_step();

	private:
		float **x;
		std::vector<std::string> take_step_block;
		std::vector<std::string> if_accept_block;
		std::vector<std::string> if_reject_block;
};
