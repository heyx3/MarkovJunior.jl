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
                @within_box inherit_biases pixel(min=1, size=3) begin
                    @rewrite (2) [R b] => [R R]
                end
            end field(R, randomness=0.3)

            # Scoped Ops run their inner Ops on a grid view, which is a different grid type.
            @within_box uv(min=0.25, size=0.5) begin
                @fill 'w' pixel(min=1, size=1)
                @rewrite (2) wb=>ww
                @sequence begin
                    @rewrite (2) [w b] => [w g]
                end field(w)
            end
            @within_regions border(0) 2 begin
                wg => @sequence begin
                    @rewrite (2) w => g
                end
                R => @fill 'G' pixel(min=1, size=1)
            end field(G)

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

    # `@grammar`, `@calc` functions, and calculation thresholds, using a condensed copy of scenes/Floorplan.jl.
    floorplan = markov_algo_parse("""
        @markovjunior 2 'L' begin
            @calc long() = grid_size_max_axis()
            @calc short() = 3 - long()
            @calc along() = size(long())
            @calc across() = size(short())
            @calc room_width() = (0.6:1.0) * across()
            @calc max_bath_width() = min(10, 0.75 * across())
            @calc n_windows() = size_average() / 4

            @grammar begin
                ROOT => inset(2, Lawn, House)
                House => inset(1, Wall, Plan)
                Plan => split(short(), WingA(w=1), Hall(3), WingB(w=1))
                WingA(along() >= across()) => repeat(long(), RoomA(each(room_width())), Wall(1))
                WingA(along() < across()) => RoomA
                WingB(along() >= across()) => repeat(long(), RoomB(each(room_width())), Wall(1))
                WingB(along() < across()) => RoomB
                RoomA(along() < max_bath_width()) => split(short(), Bath, BathWall(1))
                RoomA => split(short(), Dry, DryWall(1))
                RoomB(along() < max_bath_width()) => split(short(), BathWall(1), Bath)
                RoomB => split(short(), DryWall(1), Dry)
                DryWall => split(long(), Wall(w=(1:3)), DryDoor(2), Wall(w=(1:3)))
                BathWall => split(long(), Wall(w=(1:3)), BathDoor(1), Wall(w=(1:3)))
                Lawn => _
                Wall => 'b'
                Hall => 'w'
                Dry => 'E'
                Bath => 'T'
                DryDoor => 'Y'
                BathDoor => 'B'
            end

            @rewrite 1 [
                w b L L
                w b L L
                w b L L
            ] => [
                _ N g g
                _ N g g
                _ N g g
            ]
            @rewrite n_windows() [
                [ET] b L
                [ET] b L
            ] => [
                _ S _
                _ S _
            ]

            @within_regions border(0) begin
                EY => @sequence begin
                    @rewrite Y => N
                    @fill 'P' pixel(min=1, max=(len, 1)) +E
                    @fill 'P' pixel(min=(1, len), max=len) +E
                    @rewrite N[EP] => NY
                    @rewrite 1 begin
                        PRIORITIZE(earliest)
                        [
                            P     P    P
                            [EP] [EP] [EP]
                            [EP] [EP] [EP]
                        ] => [
                            w w w
                            R R R
                            R R R
                        ]
                        PP => NN
                    end
                    @rewrite [PY] => E
                end
            end
        end
    """)
    markov_algo_to_string(floorplan)
    channel = markov_algo_run(floorplan, (48, 28), MarkovTickSettings(), seeds=UInt32(1))
    markov_algo_complete(_ -> nothing, channel)
end