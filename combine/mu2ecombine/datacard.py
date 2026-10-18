"""A Combine datacard writer supporting several channels in one card."""


def kappa_text(down, up=None, precision=4):
    """lnN text: 'k' for a symmetric effect, 'kdown/kup' otherwise."""
    if up is None or abs(up * down - 1.) < 1e-6 and up >= 1.:
        k = up if up is not None else down
        return f"{k:.{precision}f}"
    return f"{down:.{precision}f}/{up:.{precision}f}"


class Datacard:
    def __init__(self, comment=""):
        self.comment = comment
        self.shapes = []         # (process, channel, file, nominal, systematic)
        self.bins = []           # (channel, observation)
        self.processes = []      # (channel, process, id, rate)
        self.systematics = {}    # name -> [type, {(channel, process): text}]
        self.extra = []          # lines after the systematics (rateParam, param, discrete, groups, MC statistics)

    def add_shapes(self, process, channel, file, nominal, systematic=""):
        self.shapes.append((process, channel, file, nominal, systematic))

    def add_bin(self, channel, observation):
        self.bins.append((channel, observation))

    def add_process(self, channel, process, pid, rate):
        self.processes.append((channel, process, pid, rate))

    def add_systematic(self, name, kind, channel, process, text):
        entry = self.systematics.setdefault(name, [kind, {}])
        if entry[0] != kind:
            raise ValueError(f"nuisance {name}: used as both {entry[0]} and {kind}")
        entry[1][(channel, process)] = text

    def text(self):
        lines = []
        if self.comment:
            lines += [f"# {line}" for line in self.comment.splitlines()]
        lines += ["imax *", "jmax *", "kmax *", "-" * 80]
        for s in self.shapes:
            lines.append(" ".join(x for x in ("shapes",) + s if x))
        if self.shapes:
            lines.append("-" * 80)
        w = max([len(c) for c, _ in self.bins] + [len(p) for _, p, _, _ in self.processes] + [12]) + 2
        lines.append("bin".ljust(16) + "".join(c.ljust(w) for c, _ in self.bins))
        lines.append("observation".ljust(16) + "".join(f"{o:.6g}".ljust(w) for _, o in self.bins))
        lines.append("-" * 80)
        head = 28
        lines.append("bin".ljust(head) + "".join(c.ljust(w) for c, _, _, _ in self.processes))
        lines.append("process".ljust(head) + "".join(p.ljust(w) for _, p, _, _ in self.processes))
        lines.append("process".ljust(head) + "".join(str(i).ljust(w) for _, _, i, _ in self.processes))
        lines.append("rate".ljust(head) + "".join(f"{r:.6g}".ljust(w) for _, _, _, r in self.processes))
        lines.append("-" * 80)
        for name, (kind, values) in self.systematics.items():
            row = f"{name} {kind}".ljust(head)
            row += "".join(values.get((c, p), "-").ljust(w) for c, p, _, _ in self.processes)
            lines.append(row.rstrip())
        lines += self.extra
        return "\n".join(lines) + "\n"

    def write(self, path):
        with open(path, "w") as f:
            f.write(self.text())
