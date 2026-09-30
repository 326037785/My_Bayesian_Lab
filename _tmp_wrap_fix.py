import pathlib, re
p = pathlib.Path('BayesianFilter/data_association/phd_filter.py')
s = p.read_text(encoding='utf-8-sig')
pat = re.compile(r'(?m)^(?P<ind>\s*)innovation = z - z_pred.*$')
def repl(m):
    ind = m.group('ind')
    add = '\n' + ind + 'if self.angle_wrap_idx is not None:'
    add += '\n' + ind + '    innovation[self.angle_wrap_idx] = (innovation[self.angle_wrap_idx] + np.pi) % (2 * np.pi) - np.pi'
    return m.group(0) + add
s2, n = pat.subn(repl, s)
p.write_text(s2, encoding='utf-8-sig')
print(f'wrapped {n}')
