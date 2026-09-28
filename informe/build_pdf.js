// Render informe/informe.html to informe/informe.pdf with headless Chromium (Playwright).
// Usage: node informe/build_pdf.js   (needs `npm i -g playwright` or NODE_PATH pointing to it)
const path = require('path');
const { chromium } = require('playwright');
(async () => {
  const here = __dirname;
  const browser = await chromium.launch();
  const page = await browser.newPage();
  await page.emulateMedia({ media: 'print', colorScheme: 'light' });
  await page.goto('file://' + path.join(here, 'informe.html'), { waitUntil: 'networkidle' });
  await page.pdf({
    path: path.join(here, 'informe.pdf'), format: 'A4', printBackground: true,
    displayHeaderFooter: true, headerTemplate: '<span></span>',
    footerTemplate: '<div style="width:100%;text-align:center;font:8px Arial;color:#777">' +
      '<span class="pageNumber"></span> / <span class="totalPages"></span></div>',
    margin: { top: '14mm', bottom: '14mm', left: '15mm', right: '15mm' },
  });
  await browser.close();
  console.log('wrote', path.join(here, 'informe.pdf'));
})();
