# Forward run directories

One directory per forward run, named by the user, as specified in SPEC_00 section 2.3 [D2a].

Each run directory holds its namelist at the top, the run's model ready inputs in `inputs/`,
and the products written by `casspian.forward` in `output/`. A run directory is self contained
and can be archived, copied, or diffed as a unit. An ensemble is a directory of such
directories.

Nothing is committed here yet. Run output is excluded by `.gitignore`.
