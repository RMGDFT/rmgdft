#!/usr/bin/env python3

if __name__ == '__main__':
    ref_name = 'input.110.ref.cond_12.dat'
    new_name = 'input.110.00.cond_12.dat'
  
    total_ref = 0.0
    total_new = 0.0
    col_index = 1

    with open(ref_name) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                total_ref += float(line.split()[col_index])
            except (IndexError, ValueError):
                pass

    with open(new_name) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            try:
                total_new += float(line.split()[col_index])
            except (IndexError, ValueError):
                pass
    if(abs(total_new - total_ref))/total_ref < 0.01:
        print("pass")
    else:
        print("faild")

