

class StepDescription {
	public:
		StepDescription();
		~StepDescription();
		std::vector<std::string>* get_take_step_block();

		void take_step();

	private:
		float **x;
		std::vector<std::string> take_step_block;
};
