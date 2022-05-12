import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

plt.switch_backend('agg')
plt.style.use('ggplot')

COLOR_TOP1 = '#d73027'
COLOR_TOP3 = '#fdae61'
COLOR_TOP5 = '#abd9e9'
COLOR_TOP10 = '#4575b4'

def pretty_print_category(category):
    if category == 'known_gene_known_disease': return 'Known Disease Caused by Disease-Causing Gene \n Previously Unassociated with the Disease'
    elif category == 'known_gene_disease': return 'Known Disease Caused by a Known, \n Associated Disease Gene'
    elif category == 'known_gene_new_disease': return 'Novel Disease Caused by a Gene \n Already Associated with another Disease'
    elif category == 'new_gene_known_disease': return 'Known Disease Caused by a Gene \n Previously Unassociated with any Disease'
    elif category == 'new_gene_new_disease': return 'Novel Disease Caused by a Gene \n Previously Unassociated with any Disease'
    else: return category

def calc_diff(list1, list2):
    return [l1-l2 for l1,l2 in zip(list1, list2)]

def results_to_lists(all_results_dict):
    top_1_acc, top_3_acc, top_5_acc, top_10_acc, mrr, avg_rank  = [], [], [], [], [], []
    labels = []
    for fname, results_dict in all_results_dict.items():
        labels.append(str(fname))
        top_1_acc.append(results_dict['top_1_acc'])
        top_3_acc.append(results_dict['top_3_acc'])
        top_5_acc.append(results_dict['top_5_acc'])
        top_10_acc.append(results_dict['top_10_acc'])
        mrr.append(results_dict['mrr'])
        avg_rank.append(results_dict['avg_rank'])

    return top_1_acc, top_3_acc, top_5_acc, top_10_acc, mrr, avg_rank, labels

def plot(all_results_dict, category, filename, is_udn=False, plot_labels=None):
    fig, ax = plt.subplots()
    width = 0.35

    top_1_acc, top_3_acc, top_5_acc, top_10_acc, mrr, avg_rank, labels = results_to_lists(all_results_dict)
    if plot_labels: labels = plot_labels

    ax.bar(labels, top_1_acc, width, label='Top 1', color=COLOR_TOP1) #, color='g'
    ax.bar(labels, calc_diff(top_3_acc,top_1_acc), width, bottom=top_1_acc,label='Top 3', color=COLOR_TOP3)
    ax.bar(labels, calc_diff(top_5_acc,top_3_acc), width, bottom=top_3_acc,label='Top 5', color=COLOR_TOP5)
    ax.bar(labels, calc_diff(top_10_acc, top_5_acc), width, bottom=top_5_acc,label='Top 10', color=COLOR_TOP10)
    ax.set_ylim([0,1.0])  
    ax.set_ylabel('Proportion of Patients')
    if is_udn:
        ax.set_title('UDN Patients')
    else:
        ax.set_title('Simulated Patients')

    ax.legend(title= "Causal Gene Rank")
    plt.show()
    plt.savefig(filename)

def plot_sim_vs_udn(filename, sim_results_dict, udn_results_dict):
    fig, ax = plt.subplots()
    width = 0.35
    COLOR = 'black'
    plt.rcParams['text.color'] = COLOR
    plt.rcParams['axes.labelcolor'] = COLOR
    plt.rcParams['xtick.color'] = COLOR
    plt.rcParams['ytick.color'] = COLOR

    top_1_acc, top_3_acc, top_5_acc, top_10_acc, mrr, avg_rank, labels = results_to_lists(sim_results_dict)
    top_1_acc_udn, top_3_acc_udn, top_5_acc_udn, top_10_acc_udn, mrr_udn, avg_rank_udn, labels_udn = results_to_lists(udn_results_dict)
    print(len(top_1_acc), len(top_1_acc_udn))
    x = np.arange(len(labels))
    x_udn = np.arange(len(labels_udn))

    buffer=0.07

    #color = next(ax._get_lines.prop_cycler)['color']
    top1 = ax.bar(x - width/2 - buffer, top_1_acc, width, label='Top 1', color=COLOR_TOP1) #, color='g'
    top1_udn = ax.bar(x_udn + width/2 + buffer, top_1_acc_udn, width, label='Top 1', color=COLOR_TOP1) #, color='g'

    #color = next(ax._get_lines.prop_cycler)['color']
    top3 = ax.bar(x - width/2 - buffer, calc_diff(top_3_acc,top_1_acc), width, bottom=top_1_acc,label='Top 3', color=COLOR_TOP3)
    top3_udn = ax.bar(x_udn + width/2 + buffer, calc_diff(top_3_acc_udn,top_1_acc_udn), width, bottom=top_1_acc_udn,label='Top 3', color=COLOR_TOP3)

    #color = next(ax._get_lines.prop_cycler)['color']
    top5 = ax.bar(x - width/2 - buffer, calc_diff(top_5_acc,top_3_acc), width, bottom=top_3_acc,label='Top 5', color=COLOR_TOP5)
    top5_udn = ax.bar(x_udn + width/2 + buffer, calc_diff(top_5_acc_udn,top_3_acc_udn), width, bottom=top_3_acc_udn,label='Top 5', color=COLOR_TOP5)

    #color = next(ax._get_lines.prop_cycler)['color']
    top10 = ax.bar(x - width/2 - buffer, calc_diff(top_10_acc, top_5_acc), width, bottom=top_5_acc,label='Top 10', color=COLOR_TOP10)
    top10_udn = ax.bar(x_udn + width/2 + buffer, calc_diff(top_10_acc_udn, top_5_acc_udn), width, bottom=top_5_acc_udn,label='Top 10', color=COLOR_TOP10)
    ax.set_ylim([0,1.0])  
    ax.set_xticks(x_udn)
    ax.set_xticklabels(labels_udn, y=-0.08, color='black')


    for loc in x_udn:
        ax.text(loc - width/2 - buffer, -0.05, 'Simulated', ha='center')
        ax.text(loc + width/2 + buffer, -0.05, 'UDN', ha='center')

    ax.set_ylabel('Proportion of Patients')
    ax.legend([top1, top3, top5, top10], ['Top 1', 'Top 3', 'Top 5', 'Top 10'], title= "Causal Gene Rank")

    plt.tight_layout()
    #plt.show()
    plt.savefig(filename)

def autolabel(ax, rects, mrr, top_10_acc, fontsize=8, vert_height=0.05):
    """
    Attach a text label above each bar displaying its height
    [f'{m:.2f}' for m in mrr]
    """
    for rect, m, acc in zip(rects, mrr, top_10_acc):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., acc + vert_height, #0.1 + height
                f'{m:.2f}',
                ha='center', va='bottom', fontsize=fontsize)

def plot_ablations(all_results_dict, filename, xlabel=''):
    fig, ax = plt.subplots(nrows=1, ncols=1) #, figsize=(6,8)
    width = 0.5
    FONTSIZE=9

    top_1_acc, top_3_acc, top_5_acc, top_10_acc, mrr  = [], [], [], [], []
    labels = []
    for fname, results_dict in all_results_dict.items():
        labels.append(str(fname))
        top_1_acc.append(results_dict['top_1_acc'])
        top_3_acc.append(results_dict['top_3_acc'])
        top_5_acc.append(results_dict['top_5_acc'])
        top_10_acc.append(results_dict['top_10_acc'])
        mrr.append(results_dict['mrr'])

    ax1 = ax.bar(labels, top_1_acc, width, label='Top 1', color=COLOR_TOP1) #, color='g'
    ax3 = ax.bar(labels, calc_diff(top_3_acc,top_1_acc), width, bottom=top_1_acc,label='Top 3', color=COLOR_TOP3)
    ax5 = ax.bar(labels, calc_diff(top_5_acc,top_3_acc), width, bottom=top_3_acc,label='Top 5', color=COLOR_TOP5)
    ax10 = ax.bar(labels, calc_diff(top_10_acc, top_5_acc), width, bottom=top_5_acc,label='Top 10', color=COLOR_TOP10)
    ax.set_ylim([0,1.0])  
    ax.set_ylabel('Proportion of Patients', fontsize=FONTSIZE)
    ax.set_xlabel(xlabel, fontsize=FONTSIZE)
    ax.set_xticklabels(labels, fontdict={'fontsize':6})

    autolabel(ax, ax10, mrr, top_10_acc)
    if 'ablation_phrank_results' in filename:
        leg = ax.legend(title='Diagnostic Gene Rank', prop={'size': FONTSIZE}, title_fontsize=FONTSIZE, loc='lower right') #title=


    plt.sca(ax)
    plt.xticks(rotation=90,  fontsize=FONTSIZE)
    plt.tight_layout()
    if 'ablation_phrank_results' not in filename: #-2.5
        if 'ablation_phenotypes_phrank_results' in filename: bbox_height = -1.3
        if 'gene_module_replaced' in filename: bbox_height= -1.5
        else: bbox_height = -2.5
        leg = ax.legend(bbox_to_anchor=(0.25, bbox_height), fancybox=True, title='Diagnostic Gene Rank', prop={'size': FONTSIZE}, title_fontsize=FONTSIZE) #title=

    
    #plt.show()
    plt.savefig(filename, bbox_inches='tight')
    print('plotted ablations')


def grouped_plot_all_categories(categories_fname, all_categories_results_dict, all_categories_results_dict_udn=None, sim_and_udn=False):
    fig, axes = plt.subplots(nrows=3, ncols=2, sharey=True, dpi=300, figsize=(7,6) ) #sharex=True
    axes_locs = [(0,0), (0,1), (1,0), (1,1), (2, 0), (2, 1)]

    for category, axis_ind in zip(all_categories_results_dict.keys(), axes_locs):
        if sim_and_udn:
            top1, top1_udn, top3, top3_udn, top5, top5_udn, top10, top10_udn = sim_v_udn_grouped_plot(pretty_print_category(category), fig, axes[axis_ind[0]][axis_ind[1]], all_categories_results_dict[category], all_categories_results_dict_udn[category], category, categories_fname)
        else:
            grouped_plot(pretty_print_category(category), fig, axes[axis_ind[0]][axis_ind[1]], all_categories_results_dict[category], category, categories_fname) # plot_labels
        axes[axis_ind[0]][axis_ind[1]].set_ylabel("")
        axes[axis_ind[0]][axis_ind[1]].set_title(pretty_print_category(category), fontsize = 7, pad=25)
    
    # set x axis location & rotation and remove unnecessary ones
    for ax in fig.axes:
        plt.sca(ax)
        plt.xticks(fontsize=8, rotation=0, fontweight='medium')
    axes[0][0].set_xticks([])
    axes[0][1].set_xticks([])
    axes[1][0].set_xticks([])

    # delete axis that doesn't contain a plot
    fig.delaxes(axes[2][1])

    # add a big axis, hide frame
    fig.add_subplot(111, frameon=False)
    plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
    plt.grid(b=None, which='both', axis='both')
    plt.ylabel('Proportion of Patients')

    plt.tight_layout()

    pa1 = Patch(facecolor=COLOR_TOP1, edgecolor='black')
    #pa1_udn = Patch(facecolor='#f46d43', edgecolor='black')
    pa3 = Patch(facecolor=COLOR_TOP3, edgecolor='black')
    #pa3_udn = Patch(facecolor='#fee090', edgecolor='black')
    pa5 = Patch(facecolor=COLOR_TOP5, edgecolor='black')
    #pa5_udn = Patch(facecolor='#e0f3f8', edgecolor='black')
    pa10 = Patch(facecolor=COLOR_TOP10, edgecolor='black')
    #pa10_udn = Patch(facecolor='#74add1', edgecolor='black')

    if sim_and_udn:
        #axes[2][0].legend(handles=[pa1, pa3, pa5, pa10, pa1_udn, pa3_udn, pa5_udn, pa10_udn], labels = ('', '', '', '', 'Top 1', 'Top 3', 'Top 5', 'Top 10'), ncol=2, handletextpad=0.5, handlelength=1.0, columnspacing=-0.5, bbox_to_anchor=(2.0, 0.5), fancybox=True, loc='center right', title= "Causal Gene Rank")
        axes[2][0].legend(handles=[top1, top3, top5, top10], labels = ('Top 1', 'Top 3', 'Top 5', 'Top 10'), bbox_to_anchor=(2.0, 0.5), fancybox=True, loc='center right', title= "Causal Gene Rank")
    else:
        axes[2][0].legend(bbox_to_anchor=(2.0, 0.5), fancybox=True, loc='center right', title= "Causal Gene Rank")

    plt.savefig(categories_fname)


def grouped_plot(title, fig, ax, all_results_dict, category, filename, plot_labels=None):
    width = 0.25

    top_1_acc, top_3_acc, top_5_acc, top_10_acc, mrr  = [], [], [], [], []
    labels = []
    for fname, results_dict in all_results_dict.items():
        labels.append(str(fname))
        top_1_acc.append(results_dict['top_1_acc'])
        top_3_acc.append(results_dict['top_3_acc'])
        top_5_acc.append(results_dict['top_5_acc'])
        top_10_acc.append(results_dict['top_10_acc'])
        mrr.append(results_dict['mrr'])

    if plot_labels: labels = plot_labels

    ax.bar(labels, top_1_acc, width, label='Top 1', color=COLOR_TOP1) #, color='g'
    ax.bar(labels, calc_diff(top_3_acc,top_1_acc), width, bottom=top_1_acc,label='Top 3', color=COLOR_TOP3)
    ax.bar(labels, calc_diff(top_5_acc,top_3_acc), width, bottom=top_3_acc,label='Top 5', color=COLOR_TOP5)
    ax10 = ax.bar(labels, calc_diff(top_10_acc, top_5_acc), width, bottom=top_5_acc,label='Top 10', color=COLOR_TOP10)
    autolabel(ax, ax10, mrr, top_10_acc, fontsize=6, vert_height=0.02)
    ax.set_ylim([0,1.0])  
    ax.set_ylabel('Proportion of Patients')
    ax.set_title(title, fontsize = 8)

def sim_v_udn_grouped_plot(title, fig, ax, sim_results_dict, udn_results_dict, category, filename):
    width = 0.25

    top_1_acc, top_3_acc, top_5_acc, top_10_acc, mrr, avg_rank, labels = results_to_lists(sim_results_dict)
    top_1_acc_udn, top_3_acc_udn, top_5_acc_udn, top_10_acc_udn, mrr_udn, avg_rank_udn, labels_udn = results_to_lists(udn_results_dict)
    
    x = np.arange(len(labels))
    x_udn = np.arange(len(labels_udn))

    buffer=0.07

    top1 = ax.bar(x - width/2 - buffer, top_1_acc, width, label='Top 1',  color=COLOR_TOP1) 
    top1_udn = ax.bar(x_udn + width/2 + buffer, top_1_acc_udn, width, label='Top 1', color=COLOR_TOP1) #'#f46d43'

    top3 = ax.bar(x - width/2 - buffer, calc_diff(top_3_acc,top_1_acc), width, bottom=top_1_acc,label='Top 3',  color=COLOR_TOP3)
    top3_udn = ax.bar(x_udn + width/2 + buffer, calc_diff(top_3_acc_udn,top_1_acc_udn), width, bottom=top_1_acc_udn,label='Top 3', color=COLOR_TOP3) #'#fee090'

    top5 = ax.bar(x - width/2 - buffer, calc_diff(top_5_acc,top_3_acc), width, bottom=top_3_acc,label='Top 5', color=COLOR_TOP5)
    top5_udn = ax.bar(x_udn + width/2 + buffer, calc_diff(top_5_acc_udn,top_3_acc_udn), width, bottom=top_3_acc_udn,label='Top 5', color=COLOR_TOP5) #'#e0f3f8'

    top10 = ax.bar(x - width/2 - buffer, calc_diff(top_10_acc, top_5_acc), width, bottom=top_5_acc,label='Top 10', color=COLOR_TOP10)
    top10_udn = ax.bar(x_udn + width/2 + buffer, calc_diff(top_10_acc_udn, top_5_acc_udn), width, bottom=top_5_acc_udn,label='Top 10', color=COLOR_TOP10) #'#74add1'
    
    autolabel(ax, top10, avg_rank, top_10_acc, fontsize=6, vert_height=0.02)
    autolabel(ax, top10_udn, avg_rank_udn, top_10_acc_udn, fontsize=6, vert_height=0.02)


    ax.set_ylim([0,1.0])  
    ax.set_xticks(x_udn)
    ax.set_xticklabels(labels_udn, y=-0.08, color='black')

    for loc in x_udn:
        ax.text(loc - width/2 - buffer, -0.13, 'Sim', ha='center', fontsize=6)
        ax.text(loc + width/2 + buffer, -0.13, 'UDN', ha='center', fontsize=6)

    ax.set_ylabel('Proportion of Patients')


    #ax.set_title(title, fontsize = 7, pad=20, loc='left')
    return top1, top1_udn, top3, top3_udn, top5, top5_udn, top10, top10_udn
