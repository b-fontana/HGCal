import argparse

class Argparser():
    def __init__(self):
        self.p = argparse.ArgumentParser()
        self.p.add_argument(
            '--mingenen',
            type=float,
            default=20,
            help='Minimum energy of the generated particle.'
        )
        self.p.add_argument(
            '--mingeneta',
            type=float,
            default=1.6,
            help='Minimum pseudo-rapidity (eta) of the generated particle.'
        )
        self.p.add_argument(
            '--maxgeneta',
            type=float,
            default=3.2,
            help='Maximum pseudo-rapidity (eta) of the generated particle.'
        )
        self.p.add_argument(
            '--etacuts',
            nargs='+',
            type=float,
            help='Cuts in pseudo-rapidity.'
        )
        self.p.add_argument(
            '--mask',
            type=int,
            default=0,
            help='Mask (3, 4, 5 or 6)'
        )
        self.p.add_argument(
            '--noPUFile',
            type=str,
            default='',
            help='Input file without pile-up.'
        )
        self.p.add_argument(
            '--PUFile',
            type=str,
            default='',
            help='Input file with pile-up.'
        )
        self.p.add_argument(
            '--puTag',
            type=str,
            default='',
            help='Pile-up tag.'
        )
        self.p.add_argument(
            '--plotLabel',
            type=str,
            default='Single #gamma',
            help='Pile-up tag.'
        )
        self.p.add_argument(
            '--samples',
            type=str,
            default='',
            help="Sample being analysed ('inner' or 'outer')"
        )
        self.p.add_argument(
            '--mode',
            type=int,
            default=0,
            help="Modes: 1 or 2"
        )
        self.p.add_argument(
            '--apply_weights',
            action='store_true',
            default=False,
            help='Whether to apply correction weights or not.'
        )
        self.p.add_argument(
            '--outpath',
            type=str,
            default='',
            help='Path of the directory where the outputs will be saved.'
        )
        self.FLAGS, _ = self.p.parse_known_args()

    def get_flags(self):
        return self.FLAGS

    def print_args(self):
        print("********************")
        print("Input arguments:")
        for v,k in self.FLAGS.__dict__.items():
            print("{}: {}".format(v,k))
        print("********************")
            
