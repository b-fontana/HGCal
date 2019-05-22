def add_args(parser):
    parser.add_argument(
        '--directory',
        type=str,
        default='pics/',
        help='Directory where pictures are saved.'
    )
    parser.add_argument(
        '--root_file',
        type=str,
        default='histo',
        help='Root file where the histograms were saved.'
    )
    parser.add_argument(
        '--mask',
        type=int,
        default=0,
        help='Mask applied to the geometry.'
    )
    return parser.parse_known_args()

def print_args(flags):
    print("********************")
    print("Input arguments:")
    for v,k in flags.__dict__.items():
        print("{}: {}".format(v,k))
    print("********************")
