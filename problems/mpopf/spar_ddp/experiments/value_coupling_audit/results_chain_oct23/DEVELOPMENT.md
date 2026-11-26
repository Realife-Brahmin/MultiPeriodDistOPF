# Preserved development observations

No compute job failed in this bounded continuation. Prior campaigns' failed
attempts and evidence archives remain untouched.

Before submission, synthetic local prototype tests exposed two issues. Their
raw tool outputs remain in the task transcript; the observations are recorded
here rather than presented as large10k scientific failures:

1. Lifted interleaved synthetic three-interior system reached
   COMPACT_RANK_GROWTH node=72 core=72 in the inherited cap64 factor routine.
   The new isolated chain routine therefore supports an explicit charged
   bounded local dense fallback <=192 and rejects larger rank saturation.
   No direction/value tolerance was changed. The actual large10k test used
   only local dense factors <=60, so this fallback was not needed there.
2. The synthetic randomly generated Schur boundary is not necessarily positive
   definite. Applying CG to it correctly rejected with CC_CONSUMER_NOT_SPD
   iteration=2. The unit test now tests CG on a known SPD system separately and
   explicitly verifies rejection on a negative identity. The chain matrix,
   action, vector and solve comparisons against dense remain unchanged.

Final local units and both compute-node units passed before the actual capture.
Actual direction consumer was independently validated against the unchanged
dense future before any compressed result was consumed.
