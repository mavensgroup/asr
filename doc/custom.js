// doc/custom.js
document.addEventListener("DOMContentLoaded", function() {

    var homeLink = path_to_root + "index.html";
    var asrLogoUrl = path_to_root + "media/asr.svg";

    // 1. COVER PAGE LOGIC (No Sidebar/Menu)
    if (document.getElementById("cover-page")) {
        var sidebar = document.querySelector(".sidebar");
        if (sidebar) sidebar.style.display = "none";

        var menuBar = document.querySelector(".menu-bar");
        if (menuBar) menuBar.style.display = "none";

        var content = document.querySelector(".content");
        if (content) {
            content.style.maxWidth = "100%";
            content.style.margin = "0 auto";
        }

    // 2. STANDARD PAGE LOGIC (Top asr Logo only)
    } else {
        var scrollbox = document.querySelector(".sidebar-scrollbox");

        if (scrollbox && !scrollbox.querySelector(".sidebar-logo")) {
            var link = document.createElement("a");
            link.href = homeLink;
            link.className = "sidebar-logo-link";
            link.title = "Go to Cover";

            var logo = document.createElement("img");
            logo.src = asrLogoUrl;
            logo.alt = "asr Logo";
            logo.className = "sidebar-logo";

            link.appendChild(logo);
            scrollbox.insertBefore(link, scrollbox.firstChild);
        }
    }
});

// API reference: the sidebar entry "Fortran Code Reference" (api.html) opens the FORD site (api/index.html) in a new window.
document.addEventListener("click", function(e) {
    var a = e.target.closest ? e.target.closest("a") : null;
    if (!a || !a.href) return;
    if (new URL(a.href, document.baseURI).pathname.endsWith("/api.html")) {
        e.preventDefault();
        window.open(path_to_root + "api/index.html", "_blank", "noopener");
    }
}, true);
