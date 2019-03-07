import org.gradle.internal.os.OperatingSystem
import org.jetbrains.kotlin.gradle.tasks.KotlinCompile

plugins {
    java
    kotlin("jvm") version "1.3.21"
}

group = "agustinpicard"
version = "1.0"

repositories {
    mavenCentral()
    maven { url = uri("https://dl.bintray.com/openrndr/openrndr/") }
    maven { url = uri("http://dl.bintray.com/tomasvolker/maven") }
    maven { url = uri("https://jitpack.io") }
    maven {
        url = uri("http://dl.bintray.com/kyonifer/maven")
        jcenter()
    } // koma for extra math functions
    maven { url = uri("https://mvnrepository.com/artifact/org.apache.commons/commons-math3") }
}

val openrndrVersion = "0.3.30"

val openrndrOS = when (OperatingSystem.current()) {
    OperatingSystem.WINDOWS -> "windows"
    OperatingSystem.LINUX -> "linux-x64"
    OperatingSystem.MAC_OS -> "macos"
    else -> error("unsupported OS")
}

dependencies {
    implementation(kotlin("stdlib-jdk8"))
    testCompile("junit", "junit", "4.12")

    implementation("com.github.TomasVolker:openrndr-math:d56f089253")

    implementation(group = "tomasvolker", name = "numeriko-core", version = "0.0.3")
    implementation(group = "tomasvolker", name = "kyplot", version = "0.0.1")
    implementation(group = "com.github.tomasvolker", name = "parallel-utils", version = "v1.0")

    compile(group = "com.kyonifer", name = "koma-core-ejml", version = "0.12")
    compile(group = "org.apache.commons", name = "commons-math3", version = "3.6.1")

    compile("org.openrndr:openrndr-core:$openrndrVersion")
    compile("org.openrndr:openrndr-extensions:$openrndrVersion")
    compile("org.openrndr:openrndr-ffmpeg:$openrndrVersion")

    runtime("org.openrndr:openrndr-gl3:$openrndrVersion")
    runtime("org.openrndr:openrndr-gl3-natives-$openrndrOS:$openrndrVersion")
}

configure<JavaPluginConvention> {
    sourceCompatibility = JavaVersion.VERSION_1_8
}
tasks.withType<KotlinCompile> {
    kotlinOptions.jvmTarget = "1.8"
}
