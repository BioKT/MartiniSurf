from martinisurf.pipeline import build_parser


def test_surface_mode_defaults_to_2_1():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
    ])
    assert args.surface_mode == "2-1"


def test_surface_mode_accepts_4_1():
    parser = build_parser()
    args = parser.parse_args([
        "--pdb", "1RJW",
        "--anchor", "1", "1",
        "--lx", "10",
        "--ly", "10",
        "--surface-mode", "4-1",
    ])
    assert args.surface_mode == "4-1"
