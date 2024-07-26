import os
import subprocess
import tkinter as tk
from tkinter import filedialog, messagebox

def run_scj_vote_sim(k, target_file, reference_files, output_folder):
    scj_vote_cmd = ['python', 'SCJ_Vote_Sim.py', str(k), target_file] + reference_files + [output_folder]
    result = subprocess.run(scj_vote_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # 隱藏輸出
    result.check_returncode()  # Raise an error if the command failed

def run_scj_scf(output_folder):
    scj_scf_cmd = ['python', 'SCJ_SCF.py', os.path.join(output_folder, 'tar.nondup.all'), os.path.join(output_folder, 'Best_SCJ_adjacency.txt'), output_folder]
    result = subprocess.run(scj_scf_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # 隱藏輸出
    result.check_returncode()  # Raise an error if the command failed

def run_dotplot(output_folder):
    dotplot_cmd = ['python', 'SCJ_Dotplot.py', output_folder]
    result = subprocess.run(dotplot_cmd, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)  # 隱藏輸出
    result.check_returncode()  # Raise an error if the command failed

def start_scaffolding(target_file, reference_files, output_folder):
    try:
        k = len(reference_files)  # 自動計算 k
        os.makedirs(output_folder, exist_ok=True)
        run_scj_vote_sim(k, target_file, reference_files, output_folder)
        run_scj_scf(output_folder)
        run_dotplot(output_folder)
        messagebox.showinfo("Success", f"Scaffolding completed successfully. Results are in {output_folder}")
    except Exception as e:
        messagebox.showerror("Error", str(e))

def select_file(entry):
    file_path = filedialog.askopenfilename()
    entry.delete(0, tk.END)
    entry.insert(0, file_path)

def select_files(entry):
    files = filedialog.askopenfilenames()
    entry.delete(0, tk.END)
    entry.insert(0, ' '.join(files))

def select_folder(entry):
    folder_path = filedialog.askdirectory()
    entry.delete(0, tk.END)
    entry.insert(0, folder_path)

def main():
    root = tk.Tk()
    root.title("SCJ Vote: A multiple reference-based scaffolder")
    root.geometry("700x500")  # 設置初始化視窗大小

    root.iconbitmap('icon.ico')

    # 設置窗口大小可調整時控件比例縮放
    root.grid_rowconfigure(0, weight=1)
    root.grid_rowconfigure(1, weight=1)
    root.grid_rowconfigure(2, weight=1)
    root.grid_rowconfigure(3, weight=1)
    root.grid_columnconfigure(0, weight=1)
    root.grid_columnconfigure(1, weight=3)
    root.grid_columnconfigure(2, weight=1)

    # 設置背景顏色和字型
    root.configure(bg='#f0f0f0')
    
    font_style = ('Arial', 12)
    label_font_style = ('Arial', 12, 'bold')
    
    # 設置控件
    tk.Label(root, text="Target Genome", bg='#f0f0f0', font=label_font_style).grid(row=0, column=0, sticky="ew", padx=10, pady=10)
    target_file_entry = tk.Entry(root, font=font_style)
    target_file_entry.grid(row=0, column=1, sticky="ew", padx=10, pady=10)
    tk.Button(root, text="Browse", command=lambda: select_file(target_file_entry), font=font_style, bg='#4CAF50', fg='white', relief=tk.RAISED).grid(row=0, column=2, sticky="ew", padx=10, pady=10)

    tk.Label(root, text="Reference Genomes", bg='#f0f0f0', font=label_font_style).grid(row=1, column=0, sticky="ew", padx=10, pady=10)
    reference_files_entry = tk.Entry(root, font=font_style)
    reference_files_entry.grid(row=1, column=1, sticky="ew", padx=10, pady=10)
    tk.Button(root, text="Browse", command=lambda: select_files(reference_files_entry), font=font_style, bg='#4CAF50', fg='white', relief=tk.RAISED).grid(row=1, column=2, sticky="ew", padx=10, pady=10)

    tk.Label(root, text="Output Folder", bg='#f0f0f0', font=label_font_style).grid(row=2, column=0, sticky="ew", padx=10, pady=10)
    output_folder_entry = tk.Entry(root, font=font_style)
    output_folder_entry.grid(row=2, column=1, sticky="ew", padx=10, pady=10)
    tk.Button(root, text="Browse", command=lambda: select_folder(output_folder_entry), font=font_style, bg='#4CAF50', fg='white', relief=tk.RAISED).grid(row=2, column=2, sticky="ew", padx=10, pady=10)

    tk.Button(root, text="Start Scaffolding", command=lambda: start_scaffolding(
        target_file_entry.get(),
        reference_files_entry.get().split(),
        output_folder_entry.get()
    ), font=('Arial', 14), bg='#2196F3', fg='white', relief=tk.RAISED).grid(row=3, columnspan=3, sticky="ew", padx=10, pady=20)

    root.mainloop()

if __name__ == '__main__':
    main()
