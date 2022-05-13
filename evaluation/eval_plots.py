import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch
import seaborn as sns
plt.switch_backend('agg')

COLOR_TOP1 = '#238b45'
COLOR_TOP3 = '#74c476' 
COLOR_TOP5 = '#bae4b3' 
COLOR_TOP10 =  '#ddf2d6'


#####################
# Helper Functions

def pretty_print_category(category):
    '''
    Return a longer text description of each disease-gene novelty category
    '''
    if category == 'known_gene_known_disease': return 'Known Disease Caused by Disease-Causing Gene \n Previously Unassociated with the Disease'
    elif category == 'known_gene_disease': return 'Known Disease Caused by a Known, \n Associated Disease Gene'
    elif category == 'known_gene_new_disease': return 'Novel Disease Caused by a Gene \n Already Associated with another Disease'
    elif category == 'new_gene_known_disease': return 'Known Disease Caused by a Gene \n Previously Unassociated with any Disease'
    elif category == 'new_gene_new_disease': return 'Novel Disease Caused by a Gene \n Previously Unassociated with any Disease'
    else: return category

def calc_diff(list1, list2):
    '''
    Calculate the difference between two lists
    '''
    return [l1-l2 for l1,l2 in zip(list1, list2)]

def results_to_lists(all_results_dict):
    '''
    Convert Results Dict to several lists
    '''
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

def autolabel(ax, rects, vals, top_10_acc, fontsize=8, vert_height=0.05):
    """
    Attach a text label above each bar displaying a value
    """
    for rect, m, acc in zip(rects, vals, top_10_acc):
        height = rect.get_height()
        ax.text(rect.get_x() + rect.get_width()/2., acc + vert_height, 
                f'{m:.2f}',
                ha='center', va='bottom', fontsize=fontsize)

#####################
# Plot performance for a single group of patients

def plot_top_k_acc(all_results_dict, filename, is_udn=False, plot_labels=None):
    sns.set_theme(style="whitegrid", context="talk") 

    fig, ax = plt.subplots(figsize=(12, 12)) 
    width = 0.5 

    top_1_acc, top_3_acc, top_5_acc, top_10_acc, mrr, avg_rank, labels = results_to_lists(all_results_dict)
    if plot_labels: labels = plot_labels

    print('top_1_acc', top_1_acc)
    print('labels', labels)
    ax.bar(labels, top_1_acc, width, label='Top 1', color=COLOR_TOP1)
    ax.bar(labels, calc_diff(top_3_acc,top_1_acc), width, bottom=top_1_acc,label='Top 3', color=COLOR_TOP3)
    ax.bar(labels, calc_diff(top_5_acc,top_3_acc), width, bottom=top_3_acc,label='Top 5', color=COLOR_TOP5)
    ax.bar(labels, calc_diff(top_10_acc, top_5_acc), width, bottom=top_5_acc,label='Top 10', color=COLOR_TOP10)
    ax.set_ylim([0,1.0])  
    plt.xticks(fontsize=18)
    ax.set_ylabel('Proportion of Patients', fontsize=18)
    ax.grid(axis='x')
    ax.axhline(0, color="k", clip_on=False)
    ax.legend(title= "Causal Gene Rank")


    plt.tight_layout(h_pad=1)
    sns.despine(bottom=True,left=True)
    plt.savefig(filename)
    plt.savefig(str(filename).replace('.png', '.pdf'))

#####################
# Plot performance on ablated patients

def plot_ablations(all_results_dict, filename, xlabel):
    sns.set_theme(style="whitegrid", context="talk") 

    fig, ax = plt.subplots( figsize=(12,9)) 
    width = 0.5
    FONTSIZE = 16

    top_1_acc, top_3_acc, top_5_acc, top_10_acc, mrr, avg_rank, labels = results_to_lists(all_results_dict)

    ax1 = ax.bar(labels, top_1_acc, width, label='Top 1', color=COLOR_TOP1) #, color='g'
    ax3 = ax.bar(labels, calc_diff(top_3_acc,top_1_acc), width, bottom=top_1_acc,label='Top 3', color=COLOR_TOP3)
    ax5 = ax.bar(labels, calc_diff(top_5_acc,top_3_acc), width, bottom=top_3_acc,label='Top 5', color=COLOR_TOP5)
    ax10 = ax.bar(labels, calc_diff(top_10_acc, top_5_acc), width, bottom=top_5_acc,label='Top 10', color=COLOR_TOP10)
    ax.set_ylim([0,1.0])  
    ax.set_ylabel('Proportion of Patients', fontsize=FONTSIZE)
    ax.set_xlabel(xlabel, fontsize=FONTSIZE)
    ax.set_xticklabels(labels, fontdict={'fontsize':FONTSIZE})

    autolabel(ax, ax10, avg_rank, top_10_acc, fontsize=FONTSIZE, vert_height=0.02)
    if 'ablation_phrank_results' in filename:
        leg = ax.legend(title='Diagnostic Gene Rank', prop={'size': FONTSIZE}, title_fontsize=FONTSIZE, loc='lower right') 


    plt.sca(ax)
    plt.xticks(rotation=90,  fontsize=FONTSIZE)
    plt.tight_layout()
    if 'ablation_phrank_results' not in filename: 
        if 'gene_module' in filename: bbox_height= -1.5 
        else: bbox_height = -1 
        leg = ax.legend(bbox_to_anchor=(0.25, bbox_height), fancybox=True, title='Diagnostic Gene Rank', prop={'size': FONTSIZE}, title_fontsize=FONTSIZE) 

    
    plt.savefig(filename, bbox_inches='tight')
    print('Plotted Ablations')

#####################
# Plot Performance across disease-gene novelty categories

def grouped_plot_all_categories(categories_fname, all_categories_results_dict, plot_labels=None):
    sns.set_theme(style="whitegrid", context="talk") 

    fig, axes = plt.subplots(nrows=3, ncols=2, sharey=True, dpi=300, figsize=(12,12))
    axes_locs = [(0,0), (0,1), (1,0), (1,1), (2, 0), (2, 1)]

    #print( all_categories_results_dict.keys())
    ordered_categories = [ 'known_gene_disease', 'known_gene_known_disease', 'known_gene_new_disease', 'new_gene_known_disease', 'new_gene_new_disease']
    for category, axis_ind in zip(ordered_categories, axes_locs):
        grouped_plot(pretty_print_category(category), fig, axes[axis_ind[0]][axis_ind[1]], all_categories_results_dict[category], category, categories_fname, plot_labels=plot_labels) 
        axes[axis_ind[0]][axis_ind[1]].set_ylabel("")
        axes[axis_ind[0]][axis_ind[1]].set_title(pretty_print_category(category), fontsize = 16, pad=25)
    
    # set x axis location & rotation and remove unnecessary ones
    for ax in fig.axes:
        plt.sca(ax)
        plt.xticks(fontsize=14, rotation=0, fontweight='medium')
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


    pa1 = Patch(facecolor=COLOR_TOP1, edgecolor='black')
    pa3 = Patch(facecolor=COLOR_TOP3, edgecolor='black')
    pa5 = Patch(facecolor=COLOR_TOP5, edgecolor='black')
    pa10 = Patch(facecolor=COLOR_TOP10, edgecolor='black')


    axes[2][0].legend(bbox_to_anchor=(2.0, 0.5), fancybox=True, loc='center right', title= "Causal Gene Rank")
    plt.savefig(categories_fname)

def grouped_plot(title, fig, ax, all_results_dict, category, filename, plot_labels=None):
    width = 0.25

    top_1_acc, top_3_acc, top_5_acc, top_10_acc, mrr, avg_rank, labels = results_to_lists(all_results_dict)
    if plot_labels: labels = plot_labels

    plt.tight_layout(h_pad=1)
    sns.despine(bottom=True,left=True)

    ax.bar(labels, top_1_acc, width, label='Top 1', color=COLOR_TOP1) #, color='g'
    ax.bar(labels, calc_diff(top_3_acc,top_1_acc), width, bottom=top_1_acc,label='Top 3', color=COLOR_TOP3)
    ax.bar(labels, calc_diff(top_5_acc,top_3_acc), width, bottom=top_3_acc,label='Top 5', color=COLOR_TOP5)
    ax10 = ax.bar(labels, calc_diff(top_10_acc, top_5_acc), width, bottom=top_5_acc,label='Top 10', color=COLOR_TOP10)
    autolabel(ax, ax10, avg_rank, top_10_acc, fontsize=16, vert_height=0.02)
    ax.set_ylim([0,1.0])  
    plt.xticks(fontsize=16)
    ax.grid(axis='x')
    ax.axhline(0, color="k", clip_on=False)
    ax.set_title(title, fontsize = 24)

