import pandas as pd
from pathlib import Path

# Input and output settings
excel_file = '东京攻略.xlsx'         # Your source Excel file
output_dir = Path('site')           # Folder to store generated HTML files
output_dir.mkdir(exist_ok=True)

# Write external CSS (style.css)
css_content = '''
:root {
    --base-font-size: 16px;
}
html {
    font-size: var(--base-font-size);
    /* Prevent mobile browsers from auto-adjusting text size */
    -webkit-text-size-adjust: none;
    -ms-text-size-adjust: none;
    text-size-adjust: none;
}
body {
    font-family: Arial, sans-serif;
    font-size: var(--base-font-size);
    padding: 20px;
    margin: 0;
}

/* NAVIGATION */
nav {
    background: #fff;
    border-bottom: 1px solid #eaeaea;
    font-size: var(--base-font-size);
}
nav ul {
    display: flex;
    flex-wrap: wrap;
    list-style: none;
    padding: 0;
    margin: 0;
}
nav li {
    margin: 0.5rem;
}
nav a {
    text-decoration: none;
    color: #007bff;
    padding: 0.5rem 1rem;
    border-radius: 4px;
    transition: background 0.3s;
    font-size: var(--base-font-size);
}
nav a:hover {
    background: #f0f0f0;
}

/* TABLE WRAPPER: rounding & shadow */
.table-wrapper {
    overflow-x: auto;
    margin-top: 20px;
    border-radius: 0.5rem;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1);
    font-size: var(--base-font-size);
}

/* THE TABLE */
table {
    background: #fff;
    border-collapse: separate;
    border-spacing: 0;
    /* Ensure table can expand beyond viewport if content is wide */
    width: max-content;
    font-size: var(--base-font-size);
}

/* HEADER */
thead {
    background-color: #5d5fec;
}
thead th {
    color: #fff;
    font-weight: 500;
    padding: 1rem;
    text-align: left;
    font-size: var(--base-font-size);
}

/* BODY CELLS: allow multiline, auto-height */
tbody td {
    padding: 0.75rem 1rem;
    color: #6b7280;
    white-space: pre-wrap;      /* Preserve whitespace and wrap */
    word-wrap: break-word;      /* Break long words */
    max-height: none;
    overflow: visible;
    text-overflow: clip;
    font-size: var(--base-font-size);
}

/* STRIPED ROWS */
tbody tr:nth-child(even) {
    background-color: #fafafa;
}

/* REMOVE ALL CELL BORDERS */
th, td {
    border: none;
}

/* FONT-SIZE TOGGLE BUTTON */
#font-size-button {
    position: fixed;
    top: 20px;
    right: 20px;
    background: #007bff;
    color: white;
    border: none;
    padding: 0.5rem;
    border-radius: 4px;
    cursor: pointer;
    font-size: var(--base-font-size);
}

/* Small screens: keep table wide, allow horizontal scroll */
@media (max-width: 600px) {
    body {
        padding: 10px;
    }
    #font-size-button {
        top: 10px;
        right: 10px;
        padding: 0.4rem;
    }
    /* Optionally set a minimum width larger than viewport */
    table {
        min-width: 800px;
    }
}
'''.strip()

with open(output_dir / 'style.css', 'w', encoding='utf-8') as f:
    f.write(css_content)

# JavaScript for font-size toggle with multiple options
js_snippet = '''
<script>
(function(){
    const btn = document.getElementById('font-size-button');
    const sizes = [12, 14, 16, 18, 20];  // pixel sizes for A--, A-, A, A+, A++
    const labels = ['A--', 'A-', 'A', 'A+', 'A++'];
    let index = 2; // start at 'A' (16px)
    btn.textContent = labels[index];
    btn.addEventListener('click', () => {
        index = (index + 1) % sizes.length;
        document.documentElement.style.setProperty(
            '--base-font-size',
            sizes[index] + 'px'
        );
        btn.textContent = labels[index];
    });
})();
</script>
'''.strip()

# Read all sheets and build navigation
xls = pd.ExcelFile(excel_file)
sheet_names = xls.sheet_names
nav_items = [f'<li><a href="{sheet}.html">{sheet}</a></li>' for sheet in sheet_names]
nav_html = '<nav>\n<ul>\n' + '\n'.join(nav_items) + '\n</ul>\n</nav>'

# HTML template function
def make_page(sheet_name: str, page_title: str, table_html: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">
  <title>{page_title}</title>
  <link rel="stylesheet" href="style.css">
</head>
<body>
  <button id="font-size-button">A</button>
  {nav_html}
  <h1 style="font-size: var(--base-font-size)">{page_title}</h1>
  <div class="table-wrapper">
    {table_html}
  </div>
  {js_snippet}
</body>
</html>"""

# Generate pages
for sheet in sheet_names:
    raw = pd.read_excel(excel_file, sheet_name=sheet, header=None, dtype=str)
    page_title = raw.iloc[0, 0]
    df = raw.iloc[1:].copy()
    df.columns = df.iloc[0]
    df = df.iloc[1:].reset_index(drop=True)
    df = df.applymap(lambda x: x.replace('\n', '<br>') if isinstance(x, str) else x)
    table_html = df.to_html(index=False, escape=False)
    page_content = make_page(sheet, page_title, table_html)
    (output_dir / f"{sheet}.html").write_text(page_content, encoding='utf-8')
first = sheet_names[0]
index_content = (output_dir / f"{first}.html").read_text(encoding='utf-8')
(output_dir / 'index.html').write_text(index_content, encoding='utf-8')
print(f"Generated {len(sheet_names)} pages in '{output_dir.resolve()}'")
