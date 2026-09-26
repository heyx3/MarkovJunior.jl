# Precompile a representative workload, so that the first algorithm run isn't ~13 seconds of JIT.
# Turn this off while developing, with a *LocalPreferences.toml* entry:
#     [MarkovJunior]
#     precompile_workload = false

@compile_workload begin
    # Note that we cannot load a scene from disk, as `__init__` does not run during precompilation.
    algo = markov_algo_parse("""
        @markovjunior begin
            @fill 'b'
            @rewrite 1 b=>w
            @rewrite wbb=>wgw
            @rewrite g=>w

            @fill 'b'
            @fill 'R' pixel(min=1, size=1)

            @sequence (4) begin
                @rewrite (4) [R b] => [R R]
                @rewrite (2) [b b ; R b] => [b R ; R b]
            end field(R, randomness=0.3)

            @rewrite (2) RR=>GG field(-R, combo=max)
        end
    """)
    markov_algo_to_string(algo)

    for resolution in ((16, 16), (8, 8, 8))
        for tick_min_priority in (2, 4)
            channel = markov_algo_run(algo, resolution, MarkovTickSettings(tick_min_priority), seeds=UInt32(1))
            markov_algo_complete(_ -> nothing, channel)
        end
    end
end