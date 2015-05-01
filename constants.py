"""
Constants for use in this package.
"""

# for analytics purposes, we want to test whether a population
# survived the introduction of treatment. So we check whether it
# survived more than CRASH_BUFFER cycles past the crash
CRASH_BUFFER = 25

# when simulation ends, we want to print a message
# explaining why it ended
END_POP_TOO_LARGE = "Population exceeded size limit."
END_POP_DIED_OUT = "Population died out."
END_MAX_CYCLES = "Simulation reached maximum cycle limit."

# mutation effect thresholds below which mutations will be classed as 'neutral'
# for instance if a threshold == 0.05, a mutation must change the
# relevant property by at least 5% to be considered non-neutral

# mutation rate mutations
MUT_THRESHOLD = 0.05
# beneficial mutations
BEN_THRESHOLD = 0.05
# deleterious mutations
DEL_THRESHOLD = 0.01
