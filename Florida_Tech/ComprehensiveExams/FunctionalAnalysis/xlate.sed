s/^# \(.*\)/\\item \\textbf{\1}/
s/ \*\*/ \\DEFINE{/g
s/\*\*/\}/g
s/^>/\\>/
s/^\* /    \\item /
s/^## /    \\item /
s/^### /        \\item /
s/^\[$/\\begin{align*}/
s/^\]$/\\end{align*}/
s/,d/\\,d/g
s/(\(.{-}\))/$\1$/g
s/“/``/g
s/”/"/g
