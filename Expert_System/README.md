# ES TPCM 

Exemplary implementation of an expert system for increasing energy efficiency of a throughput cleaning machine (TPCM). The system employs data-driven models to identify inefficient parameter settings, calculate achievable energy savings, and prioritize actions based on a fuzzy rule base. Proposed measures can first be applied to an analytical real-time simulation model of the production machine to verify that constraints required for the specified product quality are met. 

## Usage

The KnowledgeBase contains knowledge and facts about the TPCM. The UserInterface is the main script to identify potential energy savings. 

By executing the UserInterface in online mode (active connection to the PLC needed), the current parameter values are automatically read from the PLC , all EnPIs are calculated and finally prioritized actions are returned. Furthermore, an offline mode allows the expert system to be operated without a connection to the PLC or to test alternative parameter values manually.
