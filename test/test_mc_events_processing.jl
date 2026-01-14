using Test
using SolidStateDetectors
using Unitful
using TypedTables
using ArraysOfArrays

@testset "table_utils" begin
    # Test functions accept AbstractVector
    col = SolidStateDetectors.nested_col([[1,2], [3]], [1,2])
    @test col == [[1],[2],[3]]

    tbl = TypedTables.Table((depositions = [[1,2], [3,4]], pos = VectorOfVectors([[CartesianPoint(1u"m", 2u"m", 0u"m")], [CartesianPoint(3u"m", 4u"m", 0u"m")]])))
    split_tbl = SolidStateDetectors.split_table_by_each_charge_deposition(tbl)
    @test length(split_tbl) == 4

    translated = SolidStateDetectors.translate_event_positions(tbl, [1u"m", 2u"m", 3u"m"])
    @test all(t -> all(p -> p isa CartesianPoint, t.pos), translated)
end
