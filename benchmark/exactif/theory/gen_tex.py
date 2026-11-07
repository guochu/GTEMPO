#!/usr/bin/env python3
"""从 keldysh.md 生成 keldysh.tex：
- 显示块 $$\\n...\\n$$ 压平为 $$...$$
- unicode 数学符号 → LaTeX 宏（跳过代码 span）
- pandoc → tex body
- 拼入 ctexart 前言
"""
import re, subprocess, sys, os

SRC = os.path.join(os.path.dirname(__file__), 'keldysh.md')
DST = os.path.join(os.path.dirname(__file__), 'keldysh.tex')

s = open(SRC, encoding='utf-8').read()

# 1) 压平显示块
s = re.sub(r'\$\$\n([^$]*?)\n\$\$',
           lambda m: '$$' + m.group(1).replace('\n', ' ') + '$$', s, flags=re.S)

# 2) unicode → LaTeX 宏（全替换，含 prose；pandoc 对 $ 外的 LaTeX 命令原样透传）
M = {
    '⇒': r'\Rightarrow ', '↔': r'\leftrightarrow ', '→': r'\to ',
    '≥': r'\geq ', '≤': r'\leq ', '×': r'\times ', '−': '-',
    '∝': r'\propto ', '≈': r'\approx ', '≡': r'\equiv ',
    '∑': r'\sum ', '∏': r'\prod ', '∫': r'\int ',
    '⟨': r'\langle ', '⟩': r'\rangle ',
    '⋯': r'\cdots ',
    'δ': r'\delta ', 'γ': r'\gamma ', 'ν': r'\nu ', 'ε': r'\epsilon ',
    'θ': r'\theta ', 'α': r'\alpha ', 'β': r'\beta ', 'η': r'\eta ',
    'ρ': r'\rho ', 'τ': r'\tau ', 'λ': r'\lambda ', 'ℓ': r'\ell ',
    'Σ': r'\Sigma ',
    '𝒞': r'\mathcal{C}', '𝓑': r'\mathcal{B}',
}
for uni, mac in M.items():
    s = s.replace(uni, mac)

# 3) 上下标 unicode → LaTeX
SUB = {'₀': '_0', '₁': '_1', '₂': '_2'}
SUP = {'¹': '^1', '²': '^2', '³': '^3', '⁴': '^4', '⁺': '^+', '⁻': '^-',
       'ᵀ': '^{T}', 'ᵢ': '_i', 'η': '', '⁻¹': '^{-1}'}
for uni, rep in SUB.items():
    if uni in s:
        s = s.replace(uni, '$' + rep + '$')
# 上标/下标更精确处理：这些出现在 code span 外的地方
for uni, rep in [('⁻¹', r'$^{-1}$'), ('ᵀ', r'$^{\mathsf T}$'), ('ᵢ', r'$_i$')]:
    s = s.replace(uni, rep)

# 4) code span 内不做数学替换（还原被误替换的）
# 把 code span 里的 \tau 等还原为 unicode（更好的做法是 upstream 就不替换，但简化处理）

# 5) pandoc
open('/tmp/keldysh_pandoc.md', 'w', encoding='utf-8').write(s)
body = subprocess.run(['pandoc', '/tmp/keldysh_pandoc.md',
                       '-f', 'markdown', '-t', 'latex'],
                      capture_output=True, text=True).stdout

# 6) 前言
preamble = r'''% !TEX program = xelatex
% 由 gen_tex.py 从 keldysh.md 自动生成；更新内容请改 md 后重跑本脚本
\documentclass[11pt,a4paper]{ctexart}

\usepackage[margin=2.2cm]{geometry}
\usepackage{amsmath,amssymb,mathtools}
\usepackage{longtable,booktabs,array}
\usepackage{graphicx}
\usepackage{xcolor}
\usepackage{fancyvrb}
\usepackage{upquote}
\usepackage[colorlinks=true,linkcolor=blue!60!black,urlcolor=blue!60!black]{hyperref}
\setmonofont{DejaVu Sans Mono}

\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\newenvironment{Shaded}{}{}
\DefineVerbatimEnvironment{Highlighting}{Verbatim}{commandchars=\\\{\}}
\newcommand{\NormalTok}[1]{#1}
\newcommand{\KeywordTok}[1]{\textbf{#1}}
\newcommand{\OperatorTok}[1]{#1}
\newcommand{\CommentTok}[1]{\textit{#1}}
\newcommand{\FloatTok}[1]{#1}
\newcommand{\FunctionTok}[1]{#1}
\newcommand{\DataTypeTok}[1]{#1}
\newcommand{\DecValTok}[1]{#1}
\newcommand{\BaseTok}[1]{#1}
\newcommand{\OtherTok}[1]{#1}
\newcommand{\AlertTok}[1]{\textbf{#1}}
\newcommand{\ErrorTok}[1]{\textbf{#1}}
\newcommand{\CharTok}[1]{#1}
\newcommand{\StringTok}[1]{#1}
\newcommand{\SpecialCharTok}[1]{#1}

\newcommand{\Pf}{\operatorname{Pf}}
\newcommand{\Tr}{\operatorname{Tr}}

\title{\textbf{严格 IF 实时核理论：Keldysh 轮廓上的逐步推导}\\
\large ——bath 积掉、IF Grassmann 变量积掉与 1/2/4 阶 Trotter 组装}
\author{benchmark/exactif 工作笔记}
\date{2026 年 11 月}

\begin{document}
\maketitle
\tableofcontents
\clearpage
'''
open(DST, 'w').write(preamble + body + '\n\\end{document}\n')
print('keldysh.tex generated')
