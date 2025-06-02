import numpy as np
import matplotlib.pyplot as plt
from tabulate import tabulate # 用于格式化输出表格

# --- 精确解函数 ---
def exact_solution(x, t):
    return np.sin(np.pi * x) * np.exp(-np.pi**2 * t) + x

# --- 近似解函数 ---
def approximate_solution(x, t, coeffs_a, N):
    """
    计算近似解 theta_approx = sin(pi*x) + x + sum_{j=1}^{N} a_j(t) (x^j - x^{j+1})
    coeffs_a: 包含 a_1(t), a_2(t), ..., a_N(t) 的列表或数组
    N: 系数的个数 (即题目中的N，对应 coeffs_a 的长度)
    """
    term_sum = 0.0
    if N > 0 and len(coeffs_a) >= N: # 确保N和coeffs_a的长度一致
        for j in range(N): # j 从 0 到 N-1
            # coeffs_a[j] 对应 a_{j+1}(t)
            # (x^{j+1} - x^{j+2})
            term_sum += coeffs_a[j] * (x**(j + 1) - x**(j + 2))
    elif N == 1 and len(coeffs_a) == 1: # 特殊处理题目中N=1的情况，但函数设计上N就是系数个数
        term_sum = coeffs_a[0] * (x**1 - x**2)
    elif N == 0: # 如果传入的N是0，或者coeffs_a为空，则不加修正项
        pass
    # else:
        # 可以选择在这里打印警告或抛出异常，如果N和coeffs_a的长度不匹配
        # print(f"Warning: N ({N}) and len(coeffs_a) ({len(coeffs_a)}) mismatch in approximate_solution.")

    return np.sin(np.pi * x) + x + term_sum


# --- 均方根误差计算函数 ---
def calculate_rmse(numerical_sol, exact_sol):
    return np.sqrt(np.mean((numerical_sol - exact_sol)**2))

# --- 主程序 ---
def main():
    filepath = '../output.txt' # 数据文件路径，假设脚本在HENS_plot目录下
    time_to_evaluate = 0.20 # 所有数据都对应这个时间点

    # 读取数据
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
    except FileNotFoundError:
        print(f"错误: 文件 '{filepath}' 未找到。请确保文件路径正确。")
        print(f"当前工作目录是: {os.getcwd()}")
        print("如果脚本在 HENS_plot 目录下，output.txt 应该在项目根目录。")
        return

    data_by_method = {
        "Galerkin": [],
        "Subdomain": [],
        "Collocation": []
    }

    n_values_in_file = [1, 3, 5, 7] # 文件中每种方法对应的N值（也是coeffs_a的长度）
    method_order = ["Galerkin", "Subdomain", "Collocation"]
    current_method_idx = 0
    line_idx_in_method = 0
    temp_method_data = []

    for line_str in lines:
        if not line_str.strip():
            continue
        parts = list(map(float, line_str.strip().split()))
        time_val = parts[0]
        coeffs = parts[1:]
        # N_val_from_file 实际上是 coeffs_a 列表的长度
        # 我们也用 n_values_in_file[line_idx_in_method] 来确保对应关系
        N_val_for_coeffs = len(coeffs)
        if N_val_for_coeffs != n_values_in_file[line_idx_in_method]:
            print(f"警告: 文件中第 {sum(len(data_by_method[m]) for m in data_by_method) + len(temp_method_data) + 1} 行的系数数量 "
                  f"({N_val_for_coeffs}) 与预期的 N={n_values_in_file[line_idx_in_method]} 不匹配。将使用系数数量作为N。")
        temp_method_data.append({'time': time_val, 'coeffs': coeffs, 'N_coeffs': N_val_for_coeffs})

        line_idx_in_method += 1
        if line_idx_in_method == len(n_values_in_file):
            if current_method_idx < len(method_order):
                data_by_method[method_order[current_method_idx]] = temp_method_data
            temp_method_data = []
            line_idx_in_method = 0
            current_method_idx += 1
            if current_method_idx >= len(method_order) and line_idx_in_method == 0: #确保在所有方法处理完后停止
                break

    # 处理最后一种方法可能未完整记录的情况 (如果文件行数不是12的倍数)
    if current_method_idx < len(method_order) and temp_method_data:
         data_by_method[method_order[current_method_idx]] = temp_method_data


    # --- 绘图和误差计算 ---
    fig, axes = plt.subplots(nrows=len(n_values_in_file), ncols=len(method_order),
                             figsize=(12, 16), sharex=True, sharey=True) # 调整figsize使子图不那么扁
    fig.suptitle(f'Numerical vs. Exact Solution at t = {time_to_evaluate:.2f}', fontsize=16)

    x_points = np.linspace(0, 1, 100)
    rmse_results = [] # 用于存储所有RMSE结果以供输出

    for col_idx, method_name in enumerate(method_order):
        method_data_list = data_by_method.get(method_name, [])
        if not method_data_list:
            print(f"警告: 没有找到方法 '{method_name}' 的数据。")
            for row_idx in range(len(n_values_in_file)):
                ax = axes[row_idx, col_idx] if len(n_values_in_file) > 1 else axes[col_idx]
                ax.text(0.5, 0.5, "No Data", ha='center', va='center', fontsize=10, color='red')
                if row_idx == len(n_values_in_file) - 1 :
                     ax.set_xlabel('x')
                if col_idx == 0:
                    ax.set_ylabel(f'N={n_values_in_file[row_idx]}\nθ(x,t)', fontsize=14)
                if row_idx == 0:
                    ax.set_title(method_name)
            continue

        for row_idx, N_expected_from_file_order in enumerate(n_values_in_file):
            ax = axes[row_idx, col_idx] if len(n_values_in_file) > 1 else axes[col_idx]

            data_entry = None
            if row_idx < len(method_data_list):
                # 假设 method_data_list 中的顺序与 n_values_in_file 的顺序严格对应
                current_entry = method_data_list[row_idx]
                # 检查当前条目的系数数量是否与期望的N值匹配
                if current_entry['N_coeffs'] == N_expected_from_file_order:
                    data_entry = current_entry
                else:
                    # 如果不匹配，尝试在列表中查找第一个N_coeffs匹配的项
                    # （这可能在文件行数不完全符合预期时有用，但通常我们期望严格对应）
                    # for entry_search in method_data_list:
                    #     if entry_search['N_coeffs'] == N_expected_from_file_order:
                    #         data_entry = entry_search
                    #         break
                    print(f"警告: 方法 '{method_name}', 第 {row_idx+1} 个预期 N 值 ({N_expected_from_file_order}) "
                          f"与实际读取数据的系数数量 ({current_entry['N_coeffs']}) 不匹配。")

            if data_entry is None:
                ax.text(0.5, 0.5, f"Data Missing\nfor N={N_expected_from_file_order}", ha='center', va='center', fontsize=9, color='red')
                # 将一个标记值（如 np.nan）添加到rmse_results中，以便表格对齐
                rmse_results.append({'Method': method_name, 'N': N_expected_from_file_order, 'RMSE': np.nan})
            else:
                coeffs_a = data_entry['coeffs']
                N_for_approx = data_entry['N_coeffs'] # 使用实际的系数数量作为N

                numerical_sol_values = np.array([approximate_solution(x, time_to_evaluate, coeffs_a, N_for_approx) for x in x_points])
                exact_sol_values = exact_solution(x_points, time_to_evaluate)

                ax.plot(x_points, exact_sol_values, 'k-', label='Exact')
                ax.plot(x_points, numerical_sol_values, 'r--', label='Numerical')

                rmse = calculate_rmse(numerical_sol_values, exact_sol_values)
                rmse_results.append({'Method': method_name, 'N': N_for_approx, 'RMSE': rmse})
                ax.text(0.05, 0.05, f'RMSE: {rmse:.2e}', transform=ax.transAxes, fontsize=8,
                        bbox=dict(boxstyle='round,pad=0.3', fc='wheat', alpha=0.5))

            # 设置轴标签和标题
            if row_idx == 0:
                ax.set_title(method_name, fontsize=12) # 稍微调大标题字号
            if col_idx == 0:
                ax.set_ylabel(f'N={N_expected_from_file_order}\nθ(x,t)', fontsize=14) # 调整y轴标签字号
            if row_idx == len(n_values_in_file) - 1 :
                ax.set_xlabel('x', fontsize=10) # 调整x轴标签字号

            ax.grid(True, linestyle=':', alpha=0.7)
            if row_idx == 0 and col_idx == len(method_order) - 1:
                ax.legend(fontsize=8, loc='upper right')

    plt.tight_layout(rect=[0, 0.02, 1, 0.96]) # 调整布局，给x轴标签和主标题留空间
    plt.savefig('solution_comparison_rmse.png', dpi=300) # 提高保存图像的DPI
    plt.show()

    # --- 输出RMSE结果到控制台 ---
    print("\n--- RMSE Results ---")
    headers = ["Method", "N", "RMSE"]
    table_data = []

    # 整理数据以便tabulate输出
    # 按方法分组，然后按N排序
    output_table = {}
    for res in rmse_results:
        if res['Method'] not in output_table:
            output_table[res['Method']] = {}
        output_table[res['Method']][res['N']] = res['RMSE']

    # 构建tabulate的输入格式
    # 表头： N=1, N=3, N=5, N=7
    # 行： Method名
    header_n_values = sorted(list(n_values_in_file))
    table_headers = ["Method"] + [f"N={n_val}" for n_val in header_n_values]
    tabulate_rows = []

    for method in method_order:
        row = [method]
        if method in output_table:
            for n_val in header_n_values:
                rmse_val = output_table[method].get(n_val, np.nan) # 如果没有对应N的值，则为nan
                if np.isnan(rmse_val):
                    row.append("N/A")
                else:
                    row.append(f"{rmse_val:.3e}") # 格式化输出
        else:
            row.extend(["N/A"] * len(header_n_values)) # 如果该方法没有数据
        tabulate_rows.append(row)

    print(tabulate(tabulate_rows, headers=table_headers, tablefmt="grid"))


if __name__ == '__main__':
    import os # 引入os模块检查文件路径
    main()