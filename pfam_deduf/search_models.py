from deduf_utils import pfam_duf
import os


class pfam_model(pfam_duf):
    def __init__(self, model_dir):
        super().__init__()
        self.model_dir = model_dir
        self.model_duf = dict()

    def search_model_families(self):
        filesin = os.listdir(self.model_dir)
        count_in_model = 0
        count_not_in_model = 0
        for f in filesin:
            fname = f[:-4]
            # print(fname)
            if self.is_model_duf(fname):
                count_in_model += 1
                fpath = os.path.join(self.model_dir, f)
        count_not_in_model = len(self.list_duf) - count_in_model
        print(count_in_model, count_not_in_model)

    def is_model_duf(self, pfam):
        if pfam in self.list_duf:
            # self.model_duf[pfam] = self.list_duf[pfam]["dufid"]
            return True
        return
