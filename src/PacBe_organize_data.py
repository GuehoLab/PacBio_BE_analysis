import sys
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib import cm
import matplotlib.patches as patches
import argparse

# python3 organize_data.py lungo_output_dir target

def get_args():
    parser = argparse.ArgumentParser(description="Base Editing Substitution Analysis & Visualization")
    parser.add_argument("lungo_dir", help="Path to the LUNGO output directory")
    parser.add_argument("target", help="Target sequence to find the center position")
    return parser.parse_args()

def main():

    args = get_args()

    lungo_dir = args.lungo_dir
    target = args.target

    ref_seq = open(f'{lungo_dir}/ref_seq/ref_wo_umi.fasta').readlines()[1].upper().strip()
    cv_pos = -1

    def rc(s): return s.translate(s.maketrans('ATGC','TACG'))[::-1]

    if target in ref_seq:
        cv_pos = ref_seq.find(target)
    elif rc(target) in ref_seq:
        cv_pos = ref_seq.find(rc(target))

    print(cv_pos)

    line = open(f'{lungo_dir}/results/preprocess_count.txt').readlines()[1].split('\t')
    control_total = int(line[4])
    treated_total = int(line[10])

    cbe_subs = {}
    for n, i in enumerate(ref_seq):
        if i in 'CG':
            pos_from_cv = n - cv_pos
            if pos_from_cv >= 0 :
                pos_from_cv += 1
            if i == 'C': pat = 'CT'
            elif i == 'G': pat = 'GA'
            cbe_subs[n] = [pos_from_cv,pat,0,0,0]
    pvals = []
    window_min = 0
    window_max = 0

    with open(f'{lungo_dir}/results/mutation_pattern_p_values.txt') as f:
        next(f)
        for line in f:
            line = line.split('\t')
            mut = line[0]
            if mut.find('substitution') == -1:
                continue
            if mut.find("'C', 'T'") == -1 and mut.find("'G', 'A'") == -1:
                continue
            if mut.find("'C', 'T'") != -1:
                pat = 'CT'
            elif mut.find("'G', 'A'") != -1:
                pat = 'GA'
            pos = int(mut.split(',')[1])
            pval = float(line[4])
            pos_from_cv = pos - cv_pos
            if pos_from_cv >= 0 :
                pos_from_cv += 1
            treated_per = round(int(line[2])*100/treated_total, 2)
            control_per = round(int(line[3])*100/control_total, 2)
            cbe_subs[pos] = [pos_from_cv, pat, treated_per, control_per, pval]
            pvals.append(pval)
            if treated_per > 2:
                if pos_from_cv < window_min:
                    window_min = pos_from_cv
                if pos_from_cv > window_max:
                    window_max = pos_from_cv
    print(f'{window_min} ~ {window_max}')

    positions = list(cbe_subs.keys())
    pvals = [cbe_subs[pos][4] for pos in positions]
    reject, qvals, _, _ = multipletests(pvals, alpha=0.05, method='fdr_bh')

    for pos, qval, rej in zip(positions, qvals, reject):
        pos_from_cv, pat, treated_per, control_per, pval = cbe_subs[pos]
        cbe_subs[pos] = [pos_from_cv, pat, treated_per, control_per, pval, qval, rej]
        

    fw = open(f'{lungo_dir}/substitution_compare.txt', 'w')
    filted_res = {}
    pos_labels = []
    window_min = len(ref_seq)
    window_max = 0
    for pos in sorted(positions):
        i = cbe_subs[pos]
        fw.write(f'{pos}')
        for x in i:
            fw.write(f'\t{x}')
        fw.write('\n')
        #if i[6] == False:
        #    continue
        if i[2] > 2:
            if pos < window_min:
                window_min = pos
            if pos > window_max:
                window_max = pos
        filted_res[pos] = i
        pos_labels.append(f'{i[0]}_{i[1][0]}')

    print(f'{window_min} ~ {window_max}')
    #pos_labels   = [v[0] for v in filted_res.values()]  
    treated_per  = [v[2] for v in filted_res.values()]
    control_per  = [v[3] for v in filted_res.values()]
    qvals        = [v[5] for v in filted_res.values()]
    rejs         = [v[6] for v in filted_res.values()]


    x = range(len(pos_labels))

    width_per_x = 0.2
    fig_width = max(12, len(pos_labels) * width_per_x)  
    fig, (ax1, ax2, ax3) = plt.subplots(
        3, 1, figsize=(fig_width, 3), sharex=True
    )

    # ─────────────────────────────
    # 공통 색 범위
    vmin = min(min(treated_per), min(control_per))
    vmax = max(max(treated_per), max(control_per))

    vmin = 0
    vmax = 10

    cmap_blue = cm.get_cmap("Blues")
    cmap_red = cm.get_cmap("Reds")

    def make_heatmap(ax, values, row_label, pos_labels, vmin, vmax):
        ax.set_xticks(range(len(pos_labels)))
        ax.set_xticklabels([])
        ax.set_yticks([0.5])
        ax.set_yticklabels([row_label])

        for j, val in enumerate(values):
            cmap = cmap_blue if pos_labels[j].endswith("C") else cmap_red
            norm_val = (val - vmin) / (vmax - vmin) if vmax > vmin else 0.5
            color = cmap(norm_val)

            # 색 채우기
            ax.add_patch(plt.Rectangle((j,0), 1, 1, color=color))

            # 색 밝기 계산해서 텍스트 색 결정
            r, g, b, _ = color
            luminance = 0.299*r + 0.587*g + 0.114*b
            text_color = "white" if luminance < 0.5 else "black"
            write_str = f"{val:.2f}"
            if val < 2:
                write_str = ""

            

            # 값 쓰기
            ax.text(j+0.5, 0.5, write_str,
                    ha="center", va="center", rotation=90, 
                    fontsize=8, color=text_color, fontweight="bold")

        ax.set_xlim(0, len(pos_labels))
        ax.set_ylim(0, 1)

    # ─────────────────────────────
    # 1. Treated heatmap
    make_heatmap(ax1, treated_per, "Treated %", pos_labels, vmin, vmax)

    # 2. Control heatmap
    make_heatmap(ax2, control_per, "Control %", pos_labels, vmin, vmax)

    # 3. q-values (bar plot 그대로)
    neglog_q = [-np.log10(q) if q > 0 else 10 for q in qvals]
    ax3.bar(x, neglog_q, color="gray", align='edge')
    ax3.axhline(-np.log10(0.05), color="red", linestyle="--", label="q=0.05 cutoff")

    ax3.set_ylabel("-log(qval)")
    ax3.set_xlabel("Position")
    ax3.set_xticks([i + 0.5 for i in x])
    ax3.set_xticklabels(pos_labels, rotation=90)
    ax3.set_ylim(0, 10)
    ax3.legend()


    nums = [int(label.split('_')[0]) for label in pos_labels]

    # 강조할 region (예: 0~20)
    region_min, region_max = 0, 20

    # region에 해당하는 index 범위 구하기
    region_idx = [i for i, n in enumerate(nums) if region_min <= n <= region_max]
    if region_idx:
        xmin = min(region_idx)
        xmax = max(region_idx) + 1   # +1 해서 마지막 칸까지 포함

        rect = patches.Rectangle(
            (xmin, -0.2),           # 시작 좌표 (x, y)
            xmax - xmin,         # width
            ax1.get_ylim()[1]+0.15,    # height = y축 전체
            linewidth=3, edgecolor="dimgray", facecolor="none", 
            zorder=10            # 다른 그래프 위에 나오게
        )
        ax1.add_patch(rect)

        rect = patches.Rectangle(
            (xmin, -0.2),           # 시작 좌표 (x, y)
            xmax - xmin,         # width
            ax2.get_ylim()[1]+0.4,    # height = y축 전체
            linewidth=3, edgecolor="dimgray", facecolor="none", 
            zorder=10            # 다른 그래프 위에 나오게
        )
        ax2.add_patch(rect)

        rect = patches.Rectangle(
            (xmin, 0.15),           # 시작 좌표 (x, y)
            xmax - xmin,         # width
            ax3.get_ylim()[1]+0.2,    # height = y축 전체
            linewidth=3, edgecolor="dimgray", facecolor="none", 
            zorder=10            # 다른 그래프 위에 나오게
        )
        ax3.add_patch(rect)


    plt.tight_layout()
    plt.subplots_adjust(hspace=0)   # <-- 그래프 사이 간격 제거
    plt.savefig(lungo_dir + '/substitutions.png', bbox_inches="tight")

    print("Result saved to:" + lungo_dir + '/substitutions.png')

if __name__ == '__main__': 
    main()

        