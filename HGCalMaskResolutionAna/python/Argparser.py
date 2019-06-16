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
            default=1.4,
            help='Minimum pseudo-rapidity (eta) of the generated particle.'
        )
        self.p.add_argument(
            '--maxgeneta',
            type=float,
            default=3.2,
            help='Maximum pseudo-rapidity (eta) of the generated particle.'
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
            
