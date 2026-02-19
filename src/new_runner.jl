#############

#############

function main()
    @game_loop begin
        INIT(
            v2i(800, 800), "Markov Junior Test",
            vsync=VsyncModes.on
        )

        SETUP = begin
            algo = MarkovAlgorithm(
                CELL_CODE_BY_CHAR['b'],
                2,
                Pair{Symbol, Any}[

                ],
                AbstractMarkovOp[
                    # @rewrite 1 b=>w
                    MarkovOpRewrite1D(
                        tuple(
                            RewriteRule_Strip(
                                tuple(
                                    (CELL_CODE_BY_CHAR['b'], CELL_CODE_BY_CHAR['w'])
                                ),
                                nothing,
                                1.0f0,
                                [ GridDir(1, 1) ]
                            )
                        ),
                        1,
                        tuple(

                        )
                    ),
                    # @rewrite wbb => wgw
                    MarkovOpRewrite1D(
                        tuple(
                            RewriteRule_Strip(
                                tuple(
                                    (CELL_CODE_BY_CHAR['w'], RewriteRuleCell_Wildcard()),
                                    (CELL_CODE_BY_CHAR['b'], CELL_CODE_BY_CHAR['g']),
                                    (CELL_CODE_BY_CHAR['b'], CELL_CODE_BY_CHAR['w']),
                                ),
                                nothing,
                                1.0f0,
                                [ GridDir(1, -1), GridDir(1, 1), GridDir(2, -1), GridDir(2, 1) ]
                            )
                        ),
                        nothing,
                        tuple(

                        )
                    ),
                    # @rewrite [wg] => w
                    MarkovOpRewrite1D(
                        tuple(
                            RewriteRule_Strip(
                                tuple(
                                    (CellTypeSet('g', 'w'), CELL_CODE_BY_CHAR['w'])
                                ),
                                nothing,
                                1.0f0,
                                [ GridDir(1, 1) ]
                            )
                        ),
                        ThresholdByArea(0.25f0),
                        tuple(

                        )
                    )
                ]
            )
            algo_state = markov_algo_start(algo, (25, 25), 12345)
            algo_grid = markov_algo_grid(algo_state)

            tex = Texture(
                SimpleFormat(FormatTypes.normalized_uint, SimpleFormatComponents.RGB, SimpleFormatBitDepths.B8),
                convert(v2u, vsize(algo_grid)),
                sampler = TexSampler{2}(
                    pixel_filter=PixelFilters.rough
                )
            )
            pixel_buffer = fill(zero(v4f), size(algo_grid))
        end

        LOOP = begin
            GLFW.WindowShouldClose(LOOP.context.window) && break

            if !markov_algo_is_finished(algo, algo_state)
                markov_algo_step(algo, algo_state)
            end

            for p in one(v2i):vsize(algo_grid)
                pixel_buffer[p] = vappend(CELL_TYPES[algo_grid[p] + 1].color, @f32(1))
            end
            set_tex_color(tex, pixel_buffer)

            clear_screen(v4f(1, 0, 1, 0))
            clear_screen(@f32(1))
            simple_blit(tex)
        end

        TEARDOWN = begin
            close(tex)
            close(algo_state, algo)
        end
    end
end
